//#define EIGEN_USE_MKL_ALL
//#pragma GCC optimize("O3","unroll-loops","inline")

//#define EIGEN_USE_BLAS
//#define EIGEN_USE_MKL_VML
//#define EIGEN_USE_LAPACKE

#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <filesystem> // Require C++17
#include <regex>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <set>
#include <iterator>
#include <sys/resource.h>
#include <unistd.h>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/io/block.hpp>

//#ifdef H5_USE_EIGEN
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Easy.hpp>
//#endif

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/NumericalDiff>

#include <lbfgsb_cpp/lbfgsb.hpp>

#include <vector>
#include <set>
#include <iostream>

#include "TMatrixT.h"
#include "TH1D.h"
#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "tools.h"
#include "prob.h"
#include "ngrid.h"

#include <mfa/mfa.hpp>
#include <mfa/block_base.hpp>


//#undef basic_string_view
using namespace std;
namespace fs = std::filesystem;
using namespace std::chrono;
//using namespace lbfgsb;

#include "opts.h"

// set input and ouptut precision here, float or double
#if 0
typedef float                          real_t;
#else
typedef double                         real_t;
#endif

bool verbose=false;

 
int transformtooutergpt(int index){
      int x = int(index/(400*25));
      int y = int((index/25)%400); 
      int z = int(index%25);
      return (x*25+y)*25+z;
}

int transformtoinnergpt(int index){
      int x = int(index/(25*25))*16;
      int y = (int(index/25)%25); 
      int z = (int(index%25));
      return x*25*25+y*25+z;
}

struct FitResult {
  size_t n_iter;
  int best_grid_point;
  double best_grid_point_x;
  double best_grid_point_y;
  double best_grid_point_z;
  double last_chi_min, delta_chi;
  int flag;
  //std::vector<double> fakedataC, collspec;
};

struct FitResult2 {
  size_t n_iter;
  int best_grid_point;
  double last_chi_min, delta_chi;
  //std::vector<double> fakedataC, collspec;
};

// arguments to block foreach functions
struct DomainArgs
{
  DomainArgs(int dom_dim, int nvars) 
  {
      tot_ndom_pts = 0;
      starts.resize(dom_dim);
      ndom_pts.resize(dom_dim);
      full_dom_pts.resize(dom_dim);
      min.resize(dom_dim);
      max.resize(dom_dim);
      s.resize(nvars);
      f.resize(nvars);
      for (auto i = 0; i < nvars; i++)
      {
          s[i] = 1.0;
          f[i] = 1.0;
      }
      r = 0;
      t = 0;
      n = 0;
      multiblock = false;
      structured = true;   // Assume structured input by default
      rand_seed  = -1;
  }
  size_t              tot_ndom_pts;
  vector<int>         starts;                     // starting offsets of ndom_pts (optional, usually assumed 0)
  vector<int>         ndom_pts;                   // number of points in domain (possibly a subset of full domain)
  vector<int>         full_dom_pts;               // number of points in full domain in case a subset is taken
  vector<real_t>      min;                        // minimum corner of domain
  vector<real_t>      max;                        // maximum corner of domain
  vector<real_t>      s;                          // scaling factor for each variable or any other usage
  real_t              r;                          // x-y rotation of domain or any other usage
  vector<real_t>      f;                          // frequency multiplier for each variable or any other usage
  real_t              t;                          // waviness of domain edges or any other usage
  real_t              n;                          // noise factor [0.0 - 1.0]
  string              infile;                     // input filename
  bool                multiblock;                 // multiblock domain, get bounds from block
  bool                structured;                 // input data lies on unstructured grid
  int                 rand_seed;                  // seed for generating random data. -1: no randomization, 0: choose seed at random
};


// block
template <typename T>
struct Block : public BlockBase<T>
{
  using Base = BlockBase<T>;
  using Base::dom_dim;
  using Base::pt_dim;
  using Base::core_mins;
  using Base::core_maxs;
  using Base::bounds_mins;
  using Base::bounds_maxs;
  using Base::overlaps;
  using Base::input;
  
  static
  void* create()              { return mfa::create<Block>(); }
  
  static
  void destroy(void* b)       { mfa::destroy<Block>(b); }
  
  static
  void add(                                   // add the block to the decomposition
           int                 gid,                // block global id
           const Bounds<T>&    core,               // block bounds without any ghost added
           const Bounds<T>&    bounds,             // block bounds including any ghost region added
           const Bounds<T>&    domain,             // global data bounds
           const RCLink<T>&    link,               // neighborhood
           diy::Master&        master,             // diy master
           int                 dom_dim,            // domain dimensionality
           int                 pt_dim,             // point dimensionality
	   T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
  {
    mfa::add<Block, T>(gid, core, bounds, domain, link, master, dom_dim, pt_dim, ghost_factor);
  }
  
  static
  void save(const void* b_, diy::BinaryBuffer& bb)    { mfa::save<Block, T>(b_, bb); }
  static
  void load(void* b_, diy::BinaryBuffer& bb)          { mfa::load<Block, T>(b_, bb); }
  
  
  // read a floating point 2d scalar dataset from HDF5
  // reads masses for geometry dimension 0 from same HDF5 file
  // assigns integer values for the geometry dimension 1 from 0 to n_pts - 1
  // f = (mass, y, value)
  //template <typename V>               // type of science value being read
  void read_3d_data(
                    const       diy::Master::ProxyWithLink& cp,
                    MFAInfo&    mfa_info,
                    DomainArgs& args,
                    //Eigen::Tensor<float, 4> vals,
                    std::vector<float> vals,
                    bool  rescale)            // rescale science values
  {
    //std::cout << "$$$$ dom_dim: " << a->dom_dim << std::endl; 
    DomainArgs* a = &args;

    const int nvars       = mfa_info.nvars();
    const VectorXi mdims  = mfa_info.model_dims();

    // Resize the vectors that record error metrics
    this->max_errs.resize(nvars);
    this->sum_sq_errs.resize(nvars);
    
    // Set points per direction and compute total points in domain
    VectorXi ndom_pts(dom_dim);
    for (int i = 0; i < dom_dim; i++){
      ndom_pts(i)     =  a->ndom_pts[i];
      //std::cout << "ndom_pts: " << ndom_pts << std::endl;
    }

    // Create input data set and add to block
    //std::cout << "dom_dim, mdims, ndom_pts.prod(), ndom_pts: " << dom_dim << ", " << mdims << ", " << ndom_pts.prod() << ", " << ndom_pts << std::endl;
    input = new mfa::PointSet<T>(dom_dim, mdims, ndom_pts.prod(), ndom_pts);
    assert(vals(0) == ndom_pts(0));
    assert(vals(1) == ndom_pts(1));
    assert(vals(2) == ndom_pts(2));
    // set geometry values
    int index = 0;
    int n = 0;
    int pd = mfa_info.pt_dim();
    std::cout << "pd: " << pd << std::endl;
    for (size_t k = 0; k < (size_t)(ndom_pts(2)); k++) {
      for (size_t j = 0; j < (size_t)(ndom_pts(1)); j++) {
        for (size_t i = 0; i < (size_t)(ndom_pts(0)); i++) {
          input->domain(n, 0) = i;
          input->domain(n, 1) = j;
          input->domain(n, 2) = k;
          for (int m = 0; m < pd-dom_dim; m++) {
	    index = (i*ndom_pts(1)*ndom_pts(2) + j*ndom_pts(2) + k)*(pd-dom_dim) + m;
	    //if(i == 0 && j < 60 && k < 60 && m<1) std::cout << "index n, 3+m, i, j, k, vals: " << n << ", " << 3+m << ", " << i << ", " << j << ", " << k << ", ";
	    //if( i==0 && j==0 && m==0 ) std::cout << "index n, 3+m, i, j, k, vals: " << index << ", " << n << ", " << 3+m << ", " << i << ", " << j << ", " << k << ", ";
            input->domain(n, dom_dim+m) = vals[index];
	    //if(i == 0 && j < 60 && k < 60 && m<1) std::cout << std::setprecision(8) << vals[index] << std::endl;
	    //if( i==0 && j==0 && m==0 ) std::cout << std::setprecision(8) << vals[index] << std::endl;
            //index++;
          }	
          n++;
        }
      }
    }
    vals.clear();
    //vals.resize(0); 
    // Init params from input data (must fill input->domain first)
    input->init_params();   

    // Construct the MFA object
    this->setup_MFA(cp, mfa_info);
    
    
    // find extent of masses, values, and science variables (bins)
    // NOTE: the difference between core_mins/maxs and bounds_mins/maxs only matters
    //       if DIY is used to partition the domain space, but we set them appropriately anyway
    bounds_mins.resize(pt_dim);
    bounds_maxs.resize(pt_dim);
    core_mins.resize(dom_dim);
    core_maxs.resize(dom_dim);
    bounds_mins = input->domain.colwise().minCoeff();
    bounds_maxs = input->domain.colwise().maxCoeff();
    core_mins = bounds_mins.head(dom_dim);
    core_maxs = bounds_maxs.head(dom_dim);
    
    //std::cout << "tot_ndom_pt, input->domain(tot_ndom_pts - 1, 1), this->dom_dim = " << tot_ndom_pts << ", " << input->domain(tot_ndom_pts - 1, 1) << ", " << this->dom_dim << std::endl;
    
    // debug
    //cerr << "domain extent:\n min\n" << this->bounds_mins << "\nmax\n" << this->bounds_maxs << endl;
  }
  
};

typedef std::array<double, 3> GridPoint;

class GridPoints {
public:
  GridPoints() {}
  GridPoints(std::vector<std::vector<double>> const & m_vec_grid, int setZero=-1) {
    for (auto gp : m_vec_grid) {
      GridPoint P = {pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])};
        if (setZero>=0) P[setZero] = 0;//_gridpoints.push_back({pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});
        if(verbose) std::cout << "P[0], P[1], P[2], P[3] : " << P[0] << ", " << P[1] << ", " << P[2] << ", " << P[3] << std::endl;
        _gridpoints.push_back(P);//{pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});
    }
  }
  
  size_t NPoints()         { return _gridpoints.size(); }
  GridPoint Get(size_t index) { return _gridpoints[index]; }
  
private:
  std::vector<GridPoint> _gridpoints;
  
};

std::ostream& operator << (std::ostream& os, GridPoint const & p) {
  os << "\t0: " << p[0]  << " 1: " << p[1]  << " 2: " << p[2] <<  "\n";
  return os;
}

inline Eigen::VectorXd collapseVectorEigen(Eigen::VectorXd  const & vin, sbn::SBNconfig const & conf){
  
  // All we want is a representation with the subchannels added together
  Eigen::VectorXd cvec(conf.num_bins_total_compressed);
  cvec.setZero();
  for (int d=0; d<conf.num_detectors;++d) {
    size_t offset_in(0), offset_out(0);
    for (int i=0; i<conf.num_channels; i++) {
      size_t nbins_chan = conf.num_bins[i];
      for (int j=0; j<conf.num_subchannels[i]; j++) {
        size_t first_in   = d*conf.num_bins_detector_block            + offset_in;
        size_t first_out  = d*conf.num_bins_detector_block_compressed + offset_out;
        cvec.segment(first_out, nbins_chan).noalias() += vin.segment(first_in, nbins_chan);
        offset_in +=nbins_chan;
      }
      offset_out += nbins_chan;
    }
  }
  return cvec;
}

inline Eigen::VectorXd collapseVectorTrans(Eigen::VectorXd  const & vin, Eigen::MatrixXd  const & trans){
 
  Eigen::VectorXd cvec = trans.transpose() * vin;
  //std::cout << "full vector size, collapsed vector size: " << vin.size() << ", " << cvec.size() << std::endl;
  //std::cout << "full vector: " << vin(0) << std::endl;
  //std::cout << "coll vector: " << cvec(0) << std::endl;
  return cvec;
}

std::vector<float> getEvenElements(const std::vector<float>& vec) {
    std::vector<float> evenVector;
    int count = 0;
    for (int i = 0; i < vec.size(); ++i) {
        if (i % 2 == 0) {
            evenVector.push_back(vec[i]);
        }
    }
    return evenVector;
}

inline std::vector<float> makeSignalObject(std::string filename, std::vector<float>& dm2_index, std::vector<float>& s2_Tue_index, std::vector<float>& s2_T24_index, int & nbins ) { 

  // Open the HDF5 file for reading
  H5Easy::File file(filename.c_str(), H5Easy::File::ReadOnly);

  // Read the dataset from the HDF5 file
  // Get the spectrum
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage before spectrum: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  Eigen::MatrixXf spectrum = H5Easy::load<Eigen::MatrixXf>(file, "tree_spectrum/vec_energy_spectrum");
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after spectrum: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  // Get the dm2, s2_Tue, s2_T24
  Eigen::VectorXf dm2 = H5Easy::load<Eigen::VectorXf>(file, "tree_spectrum/value_dm2");
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after loading dm2: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  Eigen::VectorXf s2_Tue = H5Easy::load<Eigen::VectorXf>(file, "tree_spectrum/value_s2_2Tue");
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after loading s2_Tue: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  Eigen::VectorXf s2_T24 = H5Easy::load<Eigen::VectorXf>(file, "tree_spectrum/value_s2_T24");
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after loading s2_T24: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;

  // Sort the entries based on unique values
  std::set<float> dm2_unique(dm2.begin(), dm2.end());
  std::set<float> s2_Tue_unique(s2_Tue.begin(), s2_Tue.end());
  std::set<float> s2_T24_unique(s2_T24.begin(), s2_T24.end());

  // Resize dm2_index to the size of dm2_unique
  dm2_index.resize(dm2_unique.size());
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, after setting dm2_index: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  s2_Tue_index.resize(s2_Tue_unique.size());
  std::cout << "rank, after setting s2_Tue_index: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  s2_T24_index.resize(s2_T24_unique.size());
  std::cout << "rank, after setting s2_Tue_index: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  
  // Copy the elements from dm2_unique to dm2_index
  //std::copy( std::begin(dm2_unique), std::next(std::begin(dm2_unique), 20), std::begin(dm2_index));
  std::copy(dm2_unique.begin(), dm2_unique.end(), dm2_index.begin());
  std::copy(s2_Tue_unique.begin(), s2_Tue_unique.end(), s2_Tue_index.begin());
  std::copy(s2_T24_unique.begin(), s2_T24_unique.end(), s2_T24_index.begin());

  // Define grid dimensions
  //const int dim1 = dm2_unique.size(); 
  const int dim1 = dm2_index.size(); 
  const int dim2 = s2_Tue_index.size(); 
  const int dim3 = s2_T24_index.size(); 
  const int grid_size = spectrum.rows(); 
  nbins = spectrum.cols(); 
  const int total_size = dim1 * dim2 * dim3;
  std::cout << "grid_size, total_size = " << grid_size << ", " << total_size << std::endl;
  std::cout << "dim1, dim2, dim3 = " << dim1 << ", " << dim2 << ", " << dim3 << std::endl;
  // Create an Eigen tensor with size 80x60x60x1092
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, before creating tensor: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  //Eigen::Tensor<float, 4> my_grid(dim1, dim2, dim3, nbins);
  //my_grid.setZero();
  std::vector<float> my_grid(total_size*nbins);
  std::fill(my_grid.begin(), my_grid.end(), 0);

  //Initialize grid in 3d
  int gridx = -1;
  int gridy = -1;
  int gridz = -1;

  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, after creating vector: "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  // Iterate over the grid points and fill the tensor
  int index=0;
  for (int l = 0; l < grid_size; ++l) {

      auto itx  = std::find( dm2_index.begin(), dm2_index.end(), dm2(l) );
      if ( itx != dm2_index.end() ) gridx = (itx - dm2_index.begin());

      auto ity  = std::find( s2_Tue_index.begin(), s2_Tue_index.end(), s2_Tue(l) );
      if (ity != s2_Tue_index.end()) gridy = (ity - s2_Tue_index.begin());

      auto itz  = std::find( s2_T24_index.begin(), s2_T24_index.end(), s2_T24(l) );
      if (itz != s2_T24_index.end()) gridz = (itz - s2_T24_index.begin());
    
      // Find the index of the grid point in the grid vectors
      if (gridx >= 0 && gridy >= 0 && gridz >= 0) {
          for ( int bin = 0; bin < nbins; bin++ ){
	      index = (gridx*dim2*dim3 + gridy*dim3 + gridz)*nbins + bin;
	      my_grid[index] = spectrum(l,bin);
	      //if( /*( l==0 || l==30 || l==31 || l==99025 ) && bin<1 ) std::cout << "l, gridx,gridy,gridz,bin,my_grid: " <<  l << ", " << gridx << ", " << gridy << ", " << gridz << ", " << bin << ", " << std::setprecision(8) << my_grid[index] << std::endl;
	  }
      }
  }
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "index, rank, after filling tensor: "  << index << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;

  return my_grid;

}

inline void makeSignalModel(diy::mpi::communicator world, Block<real_t>* b, const diy::Master::ProxyWithLink& cp, std::vector<float> myVector, int dim1, int dim2, int dim3, int nbins, int mfadim, int deg)
{  

  // get tensor dimensions
  if(world.rank()==0) std::cout << "rank, dim1, dim2, dim3, nbins: " << world.rank() << ", " << dim1 << "," << dim2 << ", " << dim3 << ", " << nbins << std::endl;
  // default command line arguments
  int    dom_dim      = mfadim;               // dimension of domain (<= pt_dim)
  int    pt_dim       = nbins+dom_dim;        // dimension of input points
  real_t noise        = 0.0;                  // fraction of noise

  vector<int> v_degrees(dom_dim, deg);
  vector<int> v_nctrls = {dim1, dim2, dim3};

  // Info classes for construction of MFA
  ModelInfo geom_info(dom_dim);
  //std::cout << "dom_dim, pt_dim - dom_dim, v_degrees, v_nctrls: " << dom_dim << ", " << pt_dim - dom_dim << ", " << v_degrees << ", " << v_nctrls << std::endl;
  ModelInfo var_info(dom_dim, pt_dim - dom_dim, v_degrees, v_nctrls);
  MFAInfo   mfa_info(dom_dim, 1, geom_info, var_info);
  mfa_info.weighted = 0;
  //mfa_info.nvars = nbins;
  
  // echo args
  fprintf(stderr, "\n--------- Input arguments ----------\n");
  cerr <<
    "pt_dim = "         << pt_dim       << " dom_dim = "        << dom_dim      <<
        //"\ngeom_info = "  << geom_info  << " v_degrees = "    << v_degrees <<
        "\ndim1 = "  << dim1  << " dim2 = "    << dim2  << " dim3 = "    << dim3 << endl;
  
  // set input domain arguments
  DomainArgs d_args(dom_dim, pt_dim);
  d_args.n            = noise;
  d_args.min = {0, 0, 0};
  d_args.max = {dim1 - 1, dim2 - 1, dim3 - 1};
  d_args.ndom_pts[0] = dim1;
  d_args.ndom_pts[1] = dim2;
  d_args.ndom_pts[2] = dim3;

  // Get the spectrum
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage before reading : "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  double t1 = MPI_Wtime();
  double t2 = MPI_Wtime();
  b->read_3d_data(cp, mfa_info, d_args, myVector, false);
  double t3 = MPI_Wtime();
  if(world.rank()==0) std::cout << "read_3d_data took: " << t3-t1 << " seconds" << std::endl;
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after reading : "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  if(world.rank()==0) std::cout << "fixed encode block" << std::endl;
  double t4 = MPI_Wtime();
  b->fixed_encode_block(cp, mfa_info);
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after fixed encode block : "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  double t5 = MPI_Wtime();
  if(world.rank()==0) std::cout << "took: " << t5-t4 << " seconds" << std::endl;
  if(world.rank()==0) std::cout << "rangec error" << std::endl;
  double t6 = MPI_Wtime();
  b->range_error(cp, 0, true, true);
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after rangec error : "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  double t7 = MPI_Wtime();
  if(world.rank()==0) std::cout << "took: " << t7-t6 << " seconds" << std::endl;
  if(world.rank()==0) std::cout << "print block" << std::endl;
  double t8 = MPI_Wtime();
  //b->print_block(cp, 0);
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "rank, Memory usage after print block : "  << usage.ru_maxrss/1000000. << " GB" << std::endl;
  double t9 = MPI_Wtime();
  if(world.rank()==0) std::cout << "took: " << t9-t8 << " seconds" << std::endl;
  if(world.rank()==0) std::cout << "done with makeSignalModel" << std::endl;
  
}

float roundoff(float value, unsigned char prec)
{
  float pow_10 = pow(10.0f, (float)prec);
  return round(value * pow_10) / pow_10;
}


void loadData(const char* fname, std::string what, std::vector<float> & v_buffer, int & n_rows, int & n_cols) {
  H5Easy::File file(fname, H5Easy::File::ReadOnly);
  Eigen::MatrixXf _mat = H5Easy::load<Eigen::MatrixXf>(file, what);
  n_rows = _mat.rows();
  n_cols = _mat.cols();
  v_buffer = std::vector<float>(_mat.data(), _mat.data() + _mat.rows() * _mat.cols());
}

void loadData(const char* fname, std::string what, std::vector<double> & v_buffer, int & n_rows, int & n_cols) {
  H5Easy::File file(fname, H5Easy::File::ReadOnly);
  Eigen::MatrixXd _mat = H5Easy::load<Eigen::MatrixXd>(file, what);
  n_rows = _mat.rows();
  n_cols = _mat.cols();
  v_buffer = std::vector<double>(_mat.data(), _mat.data() + _mat.rows() * _mat.cols());
}

void loadData(const char* fname, std::string what, std::vector<int> & v_buffer, int & n_rows, int & n_cols) {
  H5Easy::File file(fname, H5Easy::File::ReadOnly);
  Eigen::MatrixXi _mat      = H5Easy::load<Eigen::MatrixXi>(file, what);
  n_rows = _mat.rows();
  n_cols = _mat.cols();
  v_buffer = std::vector<int>(_mat.data(), _mat.data() + _mat.rows() * _mat.cols());
}

Eigen::MatrixXd bcMatrixXd(diy::mpi::communicator world, std::vector<double>  v_buffer, int  n_rows, int  n_cols) {
  diy::mpi::broadcast(world, v_buffer, 0);
  diy::mpi::broadcast(world, n_rows,   0);
  diy::mpi::broadcast(world, n_cols,   0);
  
  Eigen::Map<Eigen::MatrixXd> mat(v_buffer.data(), n_rows, n_cols);
  return mat;
}

Eigen::MatrixXf bcMatrixXf(diy::mpi::communicator world, std::vector<float>  v_buffer, int  n_rows, int  n_cols) {
  diy::mpi::broadcast(world, v_buffer, 0);
  diy::mpi::broadcast(world, n_rows,   0);
  diy::mpi::broadcast(world, n_cols,   0);
  
  Eigen::Map<Eigen::MatrixXf> mat(v_buffer.data(), n_rows, n_cols);
  return mat;
}

Eigen::MatrixXi bcMatrixXi(diy::mpi::communicator world, std::vector<int>  v_buffer, int  n_rows, int  n_cols) {
  diy::mpi::broadcast(world, v_buffer, 0);
  diy::mpi::broadcast(world, n_rows,   0);
  diy::mpi::broadcast(world, n_cols,   0);
  
  Eigen::Map<Eigen::MatrixXi> mat(v_buffer.data(), n_rows, n_cols);
  return mat;
}

void createDataSets(HighFive::File* file, size_t nPoints, size_t nUniverses) {
  file->createDataSet<double>("last_chi_min", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<double>("delta_chi",    HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<int>("best_grid_point", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<double>("best_grid_point_x", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<double>("best_grid_point_y", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<int>("n_iter",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  // Some bookkeeping why not
  file->createDataSet<int>("i_grid",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
  file->createDataSet<int>("i_univ",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
  file->createDataSet<float>("gridx",        HighFive::DataSpace( {nPoints,        1} ));
  file->createDataSet<float>("gridy",        HighFive::DataSpace( {nPoints,        1} ));
  file->createDataSet<float>("gridz",        HighFive::DataSpace( {nPoints,        1} ));
}

void writeGrid(HighFive::File* file, std::vector<float> my_grid, std::vector<float> xcoord, std::vector<float> ycoord, std::vector<float> zcoord ) {

  HighFive::DataSet d_gridx          = file->getDataSet("gridx");
  HighFive::DataSet d_gridy          = file->getDataSet("gridy");
  HighFive::DataSet d_gridz          = file->getDataSet("gridz");
  std::cout << "x : " ;
  for(int i = 0; i < xcoord.size(); i++) std::cout << i << ", " << xcoord[i] << ", ";
  std::cout << std::endl;
  std::cout << "y : " ;
  for(int i = 0; i < ycoord.size(); i++) std::cout << i << ", " << ycoord[i] << std::endl;
  std::cout << std::endl;
  std::cout << "z : " ;
  for(int i = 0; i < zcoord.size(); i++) std::cout << i << ", " << zcoord[i] << std::endl;
  std::cout << std::endl;
  d_gridx.select(   {0, 0}, {xcoord.size(), 1}).write(xcoord);
  d_gridy.select(   {0, 0}, {ycoord.size(), 1}).write(ycoord);
  d_gridz.select(   {0, 0}, {zcoord.size(), 1}).write(zcoord);

}

Eigen::VectorXd getSpectrum(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int gpt, int dim2, int dim3 ){

    //translate to the coordinate
    Eigen::VectorXd x(3);
    x[0] = gpt / (dim2 * dim3);
    x[1] = (gpt / dim3) % dim2;
    x[2] = gpt % dim3;
    
    //prepare the MFA parameters
    int dom_dim = b->dom_dim;
    int pt_dim  = b->pt_dim;
    Eigen::VectorXd in_param(dom_dim);
    Eigen::VectorXd out_pt(pt_dim);
    Eigen::VectorXd full(pt_dim-dom_dim);
    
    //Eigen::VectorXd diff(data.size());
    // normalize the science variable to the same extent as max of geometry
    Eigen::VectorXd scale = (b->bounds_maxs - b->bounds_mins).head(b->dom_dim);
    Eigen::VectorXd scaleInv = scale.cwiseInverse();
    
    for (auto i = 0; i < dom_dim; i++){
      in_param(i) = x[i]*scaleInv(i);
    }
    //std::cout << "x, scaleInv, in_param : " << x << ", " << scaleInv << ", " << in_param << std::endl;
    
    b->decode_point(cp, in_param, out_pt);
    full = out_pt.tail(pt_dim - dom_dim);
    
    return full;
}

/*
Eigen::VectorXf getSpectrum(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int gpt, int dim2, int dim3 ){

    //translate to the coordinate
    Eigen::VectorXf x(3);
    x[0] = gpt / (dim2 * dim3);
    x[1] = (gpt / dim3) % dim2;
    x[2] = gpt % dim3;
    
    //prepare the MFA parameters
    int dom_dim = b->dom_dim;
    int pt_dim  = b->pt_dim;
    Eigen::VectorXf in_param(dom_dim);
    Eigen::VectorXf out_pt(pt_dim);
    Eigen::VectorXf full(pt_dim-dom_dim);
    
    //Eigen::VectorXf diff(data.size());
    // normalize the science variable to the same extent as max of geometry
    Eigen::VectorXf scale = (b->bounds_maxs - b->bounds_mins).head(b->dom_dim);
    Eigen::VectorXf scaleInv = scale.cwiseInverse();
    
    for (auto i = 0; i < dom_dim; i++){
      in_param(i) = x[i]*scaleInv(i);
    }
    //std::cout << "x, scaleInv, in_param : " << x << ", " << scaleInv << ", " << in_param << std::endl;
    
    b->decode_point(cp, in_param, out_pt);
    full = out_pt.tail(pt_dim - dom_dim);
    
    return full;
}
*/

inline double calcChi(Eigen::VectorXd const & data, Eigen::VectorXd const & prediction, Eigen::MatrixXd const & C_inv ) {
   auto const & diff = data-prediction;
   return diff.transpose() * C_inv * diff;
}
inline double calcChi(Eigen::VectorXd const & diff, Eigen::MatrixXd const & C_inv ) {
   //std::cout << "diff: " << diff << std::endl;
   //std::cout << "C_inv: " << C_inv.rows() << "x" <<C_inv.cols() << std::endl;
   return diff.transpose() * C_inv * diff;
}

inline double calcChi(Eigen::VectorXf const & data, Eigen::VectorXf const & prediction, Eigen::MatrixXf const & C_inv ) {
   auto const & diff = data-prediction;
   return diff.transpose() * C_inv * diff;
}
inline double calcChi(Eigen::VectorXf const & diff, Eigen::MatrixXf const & C_inv ) {
   //std::cout << "diff: " << diff << std::endl;
   //std::cout << "C_inv: " << C_inv.rows() << "x" <<C_inv.cols() << std::endl;
   return diff.transpose() * C_inv * diff;
}

double GetLLHFromSpectra(Eigen::VectorXd const & diff, Eigen::MatrixXd const & C_inv, Eigen::MatrixXd const & C_mat){
        // // function to calculate a chi2 (shape + rate)
        // // inputs:
        // // obsSpec: "data" spectra/null
        // // predSpec: "MC" spectra/gen at each point
        // // Mfracsys: total (flux+xsec+detvar) covariance matrix
        float chisqTest;
        // add the chi2-like part
	chisqTest = diff.transpose() * C_inv * diff;
        // now need ln(det(2Pi*M))
        // TMatrixD tempcov = 2*3.14159265358979323846*Msys;
	//Eigen::MatrixXd tempcov = 2*M_PI*C_mat;
	chisqTest += log(fabs(C_mat.determinant()));
	if(verbose) std::cout << "chisqTest2: " << chisqTest << std::endl;

        return chisqTest;
}//end of GetLLHFromSpectr

inline std::tuple<double, int> universeChi2(
				Block<real_t>* b, //replace with block that has encoded MFA model
				diy::Master::ProxyWithLink const& cp,
				Eigen::VectorXd const & data,
				Eigen::MatrixXd const & C_inv,
				int gridsize,
	                        int dim2,
	                        int dim3  )
{
   double chimin=std::numeric_limits<double>::infinity();
   Eigen::VectorXd diff(data.rows());
   int bestP(0);
   for (size_t i=0; i< gridsize; ++i) {
      diff = data - getSpectrum(b, cp, i, dim2, dim3);
      double chi = calcChi(diff, C_inv);
       if (chi<chimin) {
          chimin = chi;
          bestP=i;
       }
   }
   return {chimin, bestP};
}

std::vector<size_t> initialScan( Block<real_t>* b, //replace with block that has encoded MFA model
                                 diy::Master::ProxyWithLink const& cp,
		                 Eigen::VectorXd const & data, 
		                 Eigen::MatrixXd const & C_inv,
                                 double maxchi2,
	                         int gridsize,
	                         int dim2,
	                         int dim3	)
{
   std::vector<size_t> goodpoints;
   for (size_t i=0; i<gridsize; ++i) {
       double chi = calcChi(data - getSpectrum(b, cp, i, dim2, dim3), C_inv);
       if (chi<maxchi2) goodpoints.push_back(i);
   }
   return goodpoints;
}

TMatrixT<double> calcCovarianceMatrix(TMatrixT<double> const & M, std::vector<double> const & spec){
    TMatrixT<double> Mout( M.GetNcols(), M.GetNcols() );
    Mout.Zero();
    // systematics per scaled event
    for(int i =0; i<M.GetNcols(); i++) {
        for(int j =0; j<M.GetNrows(); j++) {
            Mout(i,j) = M(i,j)*spec[i]*spec[j];
            if (i==j) Mout(i,i) += spec[i];
        }
    }
    return Mout;
}

// Can we optimise this?
inline Eigen::MatrixXd calcCovarianceMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
   Eigen::MatrixXd test(M.cols(), M.cols());
    for(int i =0; i<M.cols(); i++) {
        for(int j =i; j<M.rows(); j++) {
            test(i,j) = M(i,j)*spec[i]*spec[j];
            if (i==j) test(i,i) += spec[i];
            test(j,i)=test(i,j);
        }
    }
    return test;
}
inline Eigen::MatrixXf calcCovarianceMatrix(Eigen::MatrixXf const & M, Eigen::VectorXf const & spec){
   Eigen::MatrixXf test(M.cols(), M.cols());
    for(int i =0; i<M.cols(); i++) {
        for(int j =i; j<M.rows(); j++) {
            test(i,j) = M(i,j)*spec[i]*spec[j];
            if (i==j) test(i,i) += spec[i];
            test(j,i)=test(i,j);
        }
    }
    return test;
}
inline Eigen::MatrixXd calcCovarianceMatrixFast(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec) {
   Eigen::MatrixXd ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   ret.diagonal() += spec;
   return ret;
}

inline Eigen::MatrixXf calcCovarianceMatrixFast(Eigen::MatrixXf const & M, Eigen::VectorXf const & spec) {
   Eigen::MatrixXf ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   ret.diagonal() += spec;
   return ret;
}

Eigen::MatrixXd calcMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
   Eigen::MatrixXd ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   return ret;
}

Eigen::MatrixXd calcMatrix(Eigen::MatrixXf const & M, Eigen::VectorXd const & spec){
   // Cast MatrixXf to MatrixXd
   Eigen::MatrixXd Md = M.cast<double>();
   Eigen::MatrixXd ret(Md.cols(), Md.cols());
   ret.array()    = Md.array()*(spec*spec.transpose()).array();
   return ret;
}

Eigen::MatrixXf calcMatrix(Eigen::MatrixXf const & M, Eigen::VectorXf const & spec){
   Eigen::MatrixXf ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   return ret;
}

Eigen::MatrixXd cholD(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec, double tol=1e-7) {
   auto in = calcMatrix(M, spec);
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(in);
   auto const & EV = eigensolver.eigenvalues();

   for (int i=0;i<EV.rows();++i) {
       if (EV[i]<=0) {
         if (fabs(EV[i]) < tol) for (int a=0; a<in.cols(); ++a) in(a,a) += EV[i];
       }
       if (fabs(EV[i])< tol) for (int a =0; a<in.cols(); a++) in(a,a) += tol;
   }
   Eigen::LLT<Eigen::MatrixXd> llt(in);
   return llt.matrixL();
}

Eigen::MatrixXf cholD(Eigen::MatrixXf const & M, Eigen::VectorXf const & spec, double tol=1e-7) {
   auto in = calcMatrix(M, spec);
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(in);
   auto const & EV = eigensolver.eigenvalues();

   for (int i=0;i<EV.rows();++i) {
       if (EV[i]<=0) {
         if (fabs(EV[i]) < tol) for (int a=0; a<in.cols(); ++a) in(a,a) += EV[i];
       }
       if (fabs(EV[i])< tol) for (int a =0; a<in.cols(); a++) in(a,a) += tol;
   }
   Eigen::LLT<Eigen::MatrixXf> llt(in);
   return llt.matrixL();
}

Eigen::VectorXd sample(Eigen::VectorXd const & spec, Eigen::MatrixXd const & LMAT, std::mt19937 & rng) {

   // current SBNchi uses the std::normal_distribution<float>  
   // somehow, it generated different result using std::normal_distribution<double> dist_normal(0,1) 
   // so comment out the other method for now
   std::normal_distribution<float>* m_dist_normal;
   m_dist_normal=new std::normal_distribution<float>;
   std::normal_distribution<float> dtemp(0.0,1.0);
   m_dist_normal->param(dtemp.param());
   Eigen::VectorXd RDM(spec.rows());
   //std::cout << "rng: " << rng << std::endl;
   for (int i=0;i<spec.rows();++i) RDM[i] = (*m_dist_normal)(rng);

   return LMAT*RDM + spec;
}

/*Eigen::VectorXd sample(Eigen::VectorXd const & spec, Eigen::MatrixXd const & LMAT, std::mt19937 & rng) {
   std::normal_distribution<double> dist_normal(0,1);
   Eigen::VectorXd RDM(spec.rows());
   for (int i=0;i<spec.rows();++i) RDM[i] = dist_normal(rng);
   return LMAT*RDM + spec;
}*/

Eigen::VectorXd poisson_fluctuate(Eigen::VectorXd const & spec, std::mt19937 & rng) {
   Eigen::VectorXd RDM(spec.rows());
   for (int i=0;i<spec.rows();++i) {
      std::poisson_distribution<int> dist_pois(spec[i]);
      RDM[i] = double(dist_pois(rng));
   }
   return RDM;
}

// Cholesky decomposition and solve for inverted matrix --- fastest
inline Eigen::MatrixXd invertMatrixEigen3(Eigen::MatrixXd const & M){
    return M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
}
/*
// Cholesky decomposition and solve for inverted matrix --- fastest
inline Eigen::MatrixXf invertMatrixEigen3(Eigen::MatrixXf const & M){
    return M.llt().solve(Eigen::MatrixXf::Identity(M.rows(), M.rows()));
}
*/
inline Eigen::MatrixXd collapseSubchannels(Eigen::MatrixXd const & EE, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_detector_block_compressed, conf.num_bins_detector_block_compressed);

    int mrow(0), mcol(0), mrow_out(0), mcol_out(0);
    for(int ic = 0; ic < conf.num_channels; ic++) {
        for(int jc =0; jc < conf.num_channels; jc++) {
            for(int m=0; m < conf.num_subchannels[ic]; m++) {
                for(int n=0; n< conf.num_subchannels[jc]; n++) {
                   int a, c;
                   a=mrow + n*conf.num_bins[jc];
                   c=mcol + m*conf.num_bins[ic];
                   retMat.block(mrow_out, mcol_out, conf.num_bins[jc], conf.num_bins[ic]).noalias() += EE.block(a, c, conf.num_bins[jc], conf.num_bins[ic]);
                }
            }
            mrow     += conf.num_subchannels[jc]*conf.num_bins[jc];
            mrow_out += conf.num_bins[jc];
        } // end of column loop
        mrow      = 0; // as we end this row, reSet row count, but jump down 1 column
        mrow_out  = 0;
        mcol     += conf.num_subchannels[ic]*conf.num_bins[ic];
        mcol_out += conf.num_bins[ic];
    } // end of row loop
    return retMat;
}

inline Eigen::MatrixXd collapseDetectors(Eigen::MatrixXd const & M, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_mode_block_compressed, conf.num_bins_mode_block_compressed);
    auto const & nrow = conf.num_bins_detector_block;
    auto const & crow = conf.num_bins_detector_block_compressed;
    for (int m=0; m<conf.num_detectors; m++) {
        for (int n=0; n<conf.num_detectors; n++) {
            retMat.block(n*crow, m*crow, crow, crow).noalias() = collapseSubchannels(M.block(n*nrow, m*nrow, nrow, nrow), conf);
        }
    }
    return retMat;
}

inline Eigen::MatrixXd collapseMatrix(Eigen::MatrixXd const & M, Eigen::MatrixXd const & trans){
    //Eigen::MatrixXd collapsed_matrix = trans.transpose() * M.block(0, 0, 364, 1092) * trans;
    Eigen::MatrixXd collapsed_matrix = trans.transpose() * M * trans;
    //print size of collapsed matrix
    //std::cout << "collapsed_matrix size: " << collapsed_matrix.rows() << ", " << collapsed_matrix.cols() << std::endl;
    //std::cout << "frac_matrix size: " << M.rows() << ", " << M.cols() << std::endl;
    //std::cout << "trans_matrix size: " << trans.rows() << ", " << trans.cols() << std::endl;
    return collapsed_matrix;
}
/*
inline Eigen::MatrixXf collapseMatrix(Eigen::MatrixXf const & M, Eigen::MatrixXf const & trans){
    Eigen::MatrixXf collapsed_matrix = trans.transpose() * M * trans;
    return collapsed_matrix;
}*/

Eigen::MatrixXd cholDcollapsed(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec, Eigen::MatrixXd const & trans, double tol=1e-7) {
   auto in = calcMatrix(M, spec);
   auto const & out = collapseMatrix(in, trans);
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(out);
   auto const & EV = eigensolver.eigenvalues();

   for (int i=0;i<EV.rows();++i) {
       if (EV[i]<=0) {
         if (fabs(EV[i]) < tol) for (int a=0; a<in.cols(); ++a) in(a,a) += EV[i];
       }
       if (fabs(EV[i])< tol) for (int a =0; a<in.cols(); a++) in(a,a) += tol;
   }
   Eigen::LLT<Eigen::MatrixXd> llt(out);
   return llt.matrixL();
}

inline std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> updateInvCov(Eigen::MatrixXd const & mat_frac, Eigen::MatrixXd const & mat_trans, Eigen::MatrixXd const & mat_dirt, Eigen::MatrixXd const & mat_data, Eigen::VectorXd const & spec_full) {
    auto const & cov = calcCovarianceMatrixFast(mat_frac, spec_full);
    auto const & CM = collapseMatrix(cov, mat_trans);
    auto out = CM + mat_dirt;
    return {out,invertMatrixEigen3(out)};
}


inline std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> updateInvCovCNP(Eigen::MatrixXd const & mat_frac, Eigen::MatrixXd const & mat_trans, Eigen::MatrixXd const & mat_dirt, Eigen::MatrixXd const & mat_data, Eigen::VectorXd const & spec_full, Eigen::VectorXd const & data, std::ofstream& debugFile) {

    debugFile << "fractional matrix (syst only): [" << mat_frac.diagonal().transpose() << "]" << std::endl;
    auto const & cov = calcMatrix(mat_frac, spec_full);
    debugFile << "full matrix (syst only): [" << cov.diagonal().transpose() << "]" << std::endl;
    auto out = collapseMatrix(cov, mat_trans);//collapse before adding on CNP terms
    debugFile << "collapsed matrix (syst only): [" << out.diagonal().transpose() << "]" << std::endl;
    out += mat_dirt;
    debugFile << "collapsed matrix (syst + dirt): [" << out.diagonal().transpose() << "]" << std::endl;
    Eigen::VectorXd const & cnp = 3.0*( data.array().inverse()+2.0*(collapseVectorTrans(spec_full, mat_trans)).array().inverse()).inverse();
    out.diagonal() += cnp;
    debugFile << "collapsed matrix (syst + dirt + CNP stats): [" << out.diagonal().transpose() << "]" << std::endl;
    return {out,invertMatrixEigen3(out)};
}
inline std::tuple<Eigen::MatrixXd,Eigen::MatrixXd> updateInvCovCNP(Eigen::MatrixXd const & mat_frac, Eigen::MatrixXd const & mat_trans, Eigen::MatrixXd const & mat_dirt, Eigen::MatrixXd const & mat_data, Eigen::VectorXd const & spec_full, Eigen::VectorXd const & data) {

    auto const & cov = calcMatrix(mat_frac, spec_full);
    auto out = collapseMatrix(cov, mat_trans);//collapse before adding on CNP terms
    out += mat_dirt;
    Eigen::VectorXd const & cnp = 3.0*( data.array().inverse()+2.0*(collapseVectorTrans(spec_full, mat_trans)).array().inverse()).inverse();
    out.diagonal() += cnp;
    return {out,invertMatrixEigen3(out)};
}
/*
inline std::tuple<Eigen::MatrixXf,Eigen::MatrixXf> updateInvCovCNP(Eigen::MatrixXf const & mat_frac, Eigen::MatrixXf const & mat_trans, Eigen::MatrixXf const & mat_dirt, Eigen::MatrixXf const & mat_data, Eigen::VectorXd const & spec_full, Eigen::VectorXd const & data) {
    auto const & cov = calcMatrix(mat_frac, spec_full);
    auto out = collapseMatrix(cov, mat_trans);//collapse before adding on CNP terms
    out += mat_dirt;
    Eigen::VectorXd const & cnp = 3.0*( data.array().inverse()+2.0*(collapseVectorTrans(spec_full, mat_trans)).array().inverse()).inverse();
    out.diagonal() += cnp;
    return {out,invertMatrixEigen3(out)};
}

inline std::tuple<Eigen::MatrixXf,Eigen::MatrixXf> updateInvCov(Eigen::MatrixXf const & mat_frac, Eigen::MatrixXf const & mat_trans, Eigen::MatrixXf const & mat_dirt, Eigen::MatrixXf const & mat_data, Eigen::VectorXf const & spec_full) {
    auto const & cov = calcCovarianceMatrixFast(mat_frac, spec_full);
    auto const & CM = collapseMatrix(cov, mat_trans);
    auto out = CM + mat_dirt;
    return {out,invertMatrixEigen3(out)};
}
*/
inline std::tuple<double, int, double> universeChi2( Block<real_t>* b, //replace with block that has encoded MFA model
                                                     diy::Master::ProxyWithLink const& cp, 
						     Eigen::VectorXd const & data, 
						     Eigen::MatrixXd const & covmat, 
						     Eigen::MatrixXd const & mat_trans, 
						     Eigen::MatrixXd const & mat_dirt, 
						     Eigen::MatrixXd const & mat_data, 
						     size_t const & in_grid,
	                                             int gridsize,
	                                             int dim2,
	                                             int dim3	)
{

    double chimin=std::numeric_limits<double>::infinity();
    Eigen::VectorXd diff(data.rows());
    int bestP(0);
    double this_chi = -999.;
    double chi = -999.;
    for (size_t i=0; i<gridsize; ++i) {
        //Calculate current full covariance matrix, collapse it, then Invert.
        auto const & temp_spec  = getSpectrum(b, cp, i, dim2, dim3);
        auto const & rescnp = updateInvCovCNP(covmat, mat_trans, mat_dirt, mat_data, temp_spec, data);
        auto const & invcov = std::get<1>(rescnp);
        auto const & covcol = std::get<0>(rescnp);
        diff = data - collapseVectorTrans(temp_spec,mat_trans);
        chi = calcChi(diff, invcov);
        if (chi<chimin) {
            chimin = chi;
            bestP=i;
        }
        if(i==in_grid){
          this_chi=chi; 
        }
        
    }
    return {chimin, bestP,this_chi};
}

inline std::vector<Eigen::MatrixXd> calcTotalMatrixDerivatives( Eigen::MatrixXd const & mat_frac, Eigen::MatrixXd const & mat_trans, Eigen::MatrixXd const & mat_dirt, Eigen::MatrixXd const & mat_data, Eigen::VectorXd const &  data, Eigen::VectorXd const & coll, int dom_dim, Eigen::VectorXd const & scaleInv, Eigen::VectorXd const & full, Eigen::MatrixXd const & fullderiv, Eigen::MatrixXd coll_deriv, bool pearson ){

      //compute derivative of the new chi2 definition by its parts
      //Mtot = Msys + Mstat
      
      //derivative of Msys: S(q) * F_ij * S(q) =  d_k{S(q)} * F_ij * S(q) + S(q) * F_ij * d_k{S(q)}
      //Here, CM: F_ij, full: full vector of S(q), coll: collapsed vector of S(q), data: fluctuated spectrum-->fixed
      //out_pt_deriv_mat: gradient of S(q), scaleInv: scale to map (u,v) space back to (x,y)
      Eigen::MatrixXd dfullderivmat(full.size(),full.size());
      Eigen::MatrixXd dfullmat(full.size(),full.size());
      Eigen::MatrixXd dMstat(data.size(),data.size());
      Eigen::MatrixXd dMsys(full.size(),full.size());

      //vector to store the derivative of the matrix in each dimension
      std::vector<Eigen::MatrixXd> dMtot;

      for (int i = 0; i < dom_dim; i++)
      {

	 //fill the diagonal of a matrix with the d_k{S(q)}, do the same for S(q)
	 dfullderivmat.setZero();
	 dfullderivmat.diagonal() += fullderiv.col(i);
	 dfullmat.setZero();
	 dfullmat.diagonal() += full;

         auto const & dMsysfull = 2*(dfullderivmat * mat_frac * dfullmat);

         //calculate the derivative of the CNP vector
	 Eigen::VectorXd const & dcnp = 6.0*data.array()*data.array()*(( 2.0*data.array() + coll.array())*(2.0*data.array() + coll.array())).inverse()*coll_deriv.col(i).array();
	 //initialize and fill diagonal of stat error matrix with the derivative of the CNP term
	 dMstat.setZero();
	 if ( !pearson ) dMstat.diagonal() += dcnp;

	 //push the matrix to vector
	 dMtot.push_back(dMsys + dMstat);
      }

      return dMtot;
}

inline std::vector<Eigen::MatrixXf> calcTotalMatrixDerivatives( Eigen::MatrixXf const & mat_frac, Eigen::MatrixXf const & mat_trans, Eigen::MatrixXf const & mat_dirt, Eigen::MatrixXf const & mat_data, Eigen::VectorXf const &  data, Eigen::VectorXf const & coll, int dom_dim, Eigen::VectorXf const & scaleInv, Eigen::VectorXf const & full, Eigen::MatrixXf const & fullderiv, Eigen::MatrixXf coll_deriv, bool pearson ){

      //compute derivative of the new chi2 definition by its parts
      //Mtot = Msys + Mstat
      
      //derivative of Msys: S(q) * F_ij * S(q) =  d_k{S(q)} * F_ij * S(q) + S(q) * F_ij * d_k{S(q)}
      //Here, CM: F_ij, full: full vector of S(q), coll: collapsed vector of S(q), data: fluctuated spectrum-->fixed
      //out_pt_deriv_mat: gradient of S(q), scaleInv: scale to map (u,v) space back to (x,y)
      Eigen::MatrixXf dfullderivmat(full.size(),full.size());
      Eigen::MatrixXf dfullmat(full.size(),full.size());
      Eigen::MatrixXf dMstat(data.size(),data.size());
      Eigen::MatrixXf dMsys(full.size(),full.size());

      //vector to store the derivative of the matrix in each dimension
      std::vector<Eigen::MatrixXf> dMtot;

      for (int i = 0; i < dom_dim; i++)
      {

	 //fill the diagonal of a matrix with the d_k{S(q)}, do the same for S(q)
	 dfullderivmat.setZero();
	 dfullderivmat.diagonal() += fullderiv.col(i);
	 dfullmat.setZero();
	 dfullmat.diagonal() += full;

         auto const & dMsysfull = 2*(dfullderivmat * mat_frac * dfullmat);

         //calculate the derivative of the CNP vector
	 Eigen::VectorXf const & dcnp = 6.0*data.array()*data.array()*(( 2.0*data.array() + coll.array())*(2.0*data.array() + coll.array())).inverse()*coll_deriv.col(i).array();
	 //initialize and fill diagonal of stat error matrix with the derivative of the CNP term
	 dMstat.setZero();
	 if ( !pearson ) dMstat.diagonal() += dcnp;

	 //push the matrix to vector
	 dMtot.push_back(dMsys + dMstat);
      }

      return dMtot;
}
/*
namespace lbfgsb{
  class LLR{
  private:
    Block<real_t>&     b;            //block encoded with MFA model
    diy::Master::ProxyWithLink const& cp;
    VectorXf    data;
    MatrixXf    mat_frac;
    MatrixXf    mat_trans;
    MatrixXf    mat_dirt;
    MatrixXf    mat_data;
    int         dim1;
    int         dim2;
    int         dim3;
    bool        pearson;
    int         dom_dim;
    int         pt_dim;
    int         dim;
    VectorXf    param;
    VectorXf    out_pt;
    VectorXf    out_pt_deriv;
    VectorXf    full;
    VectorXf    coll;
    VectorXf    diff;
    VectorXf    grad_eigen;
    MatrixXf    out_pt_deriv_mat;
    MatrixXf    coll_deriv;
    MatrixXf    full_deriv;
    VectorXf    scale;        // vector containing size of domain in each dimension
    VectorXf    scaleInv;     // elementwise reciprocal of scale (to cut down on FP divisions)
    int         itersgrad{0};
  public:
    //need mfa model, fake_data, and inverse covariance matrix
    LLR( Block<real_t>& b_,  
    diy::Master::ProxyWithLink const& cp_,
    VectorXf data_,
    MatrixXf mat_frac_,
    MatrixXf mat_trans_,
    MatrixXf mat_dirt_,
    MatrixXf mat_data_,
    int dim1_,
    int dim2_,
    int dim3_,
    bool pearson_,
    int dim_   ) : 
        b(b_), 	
        cp(cp_),
        data(data_),
        mat_frac(mat_frac_),
        mat_trans(mat_trans_),
        mat_dirt(mat_dirt_),
        mat_data(mat_data_),
        dim1(dim1_),
        dim2(dim2_),
        dim3(dim3_),
        pearson(pearson_),
        dom_dim(b.dom_dim),
        pt_dim(b.pt_dim),
        dim(dim_),
        param(dom_dim),
        out_pt(pt_dim),
        out_pt_deriv(pt_dim),
        full(pt_dim-dom_dim),
        coll(data.size()),
        diff(data.size()),
        grad_eigen(dom_dim),
        out_pt_deriv_mat(pt_dim-dom_dim, dom_dim),
        coll_deriv(data.size(), dom_dim),
        full_deriv(pt_dim-dom_dim, dom_dim)
    {
      scale = (b.bounds_maxs - b.bounds_mins).head(b.dom_dim);
      scaleInv = scale.cwiseInverse();
    }
    //TO DO: Figure out the correct way to pass the MFA 
    //using typename BoundedProblem<T>::TVector;

    double operator()(std::array<double, 3>& x, std::array<double, 3>& grad){
      itersgrad++;
      std::vector<double> bds; bds.push_back(dim1-1); bds.push_back(dim2-1); bds.push_back(dim3-1);
      for (int i = 0; i < dom_dim; i++){
        if (x[i] < 0 && x[i])
        {
          //std::cout << /*std::setprecision(18) <<* /  "Point value less than 0!" << std::endl;
          //std::cout << /*std::setprecision(18) <<* / "x[" << i << "] = " << x[i] << std::endl;
	  x[i] = 0.;
          //exit(1);
        }
        if (x[i] > bds[i])
        {
          //std::cout << /*std::setprecision(18) <<* / "Point value greater than upper bound!" << std::endl;
          //std::cout << /*std::setprecision(18) <<* / "x[" << i << "] = " << x[i] << std::endl;
	  x[i] = bds[i];
          //exit(1);
        }
      }

      for (auto i = 0; i < dom_dim; i++){
        param(i) = x[i]*scaleInv(i);
      }

      for (int i = 0; i < dom_dim; i++)
      {
        if (param(i) < 0)
        {
          //std::cout << /*std::setprecision(18) <<* / "Param value less than 0!" << std::endl;
          //std::cout << /*std::setprecision(18) <<* / "u[" << i << "] = " << param(i) << std::endl;
	  param(i) = 0.;
          //exit(1);
        }
        if (param(i) > 1)
        {
          //std::cout << /*std::setprecision(18) <<* / "Param value greater than 1!" << std::endl;
          //std::cout << /*std::setprecision(18) <<* / "u[" << i << "] = " << param(i) << std::endl;
	  param(i) = 1.0;
          //exit(1);
        }
      }

      //decode point from the MFA model
      b.decode_point(cp, param, out_pt);
      full = out_pt.tail(pt_dim - dom_dim);
      
      //collapse:
      coll = collapseVectorTrans(full, mat_trans);   // necessary? full.size() = nBins = pt_dim-dom_dim already

      Eigen::MatrixXf invcov;
      Eigen::MatrixXf covcol;

      if( !pearson ) {
        //get covariance matrix and inverse covariance matrix after adding CNP stats
        auto const & rescnp = updateInvCovCNP(mat_frac, mat_trans, mat_dirt, mat_data, full, data);
        invcov   = std::get<1>(rescnp);
        covcol   = std::get<0>(rescnp);
      } else {
        //get covariance matrix and inverse covariance matrix with pearson format
        auto const & rescnp = updateInvCov(mat_frac, mat_trans, mat_dirt, mat_data, full);
        invcov   = std::get<1>(rescnp);
        covcol   = std::get<0>(rescnp);
      }

      // compute objective function
      diff = data - coll;
      double fun = calcChi(diff, invcov);

      //find gradient
      for (auto i = 0; i < dom_dim; i++){
	      b.differentiate_point(cp, param, 1, i, -1, out_pt_deriv);
	      out_pt_deriv_mat.col(i) = out_pt_deriv.tail(pt_dim - dom_dim);
      }

      // scale the derivatives from (u,v) space back to x,y
      // collapse
      for (auto q = 0; q < dom_dim; q++){
        full_deriv.col(q) = 2*scaleInv(q) * out_pt_deriv_mat.col(q);
        coll_deriv.col(q) = collapseVectorTrans(full_deriv.col(q), mat_trans); // necessary?
      }
      auto const & dMtot = calcTotalMatrixDerivatives( mat_frac, mat_trans, mat_dirt, mat_data, data, coll, dom_dim, scaleInv, full, full_deriv, coll_deriv, pearson);

      // calculate the new function
      // compute gradient of objective function
      auto const & v = invcov*diff;
      for (auto i = 0; i < dom_dim; i++){
         grad_eigen(i) = -1.0*v.transpose()*dMtot[i]*v;
         grad_eigen(i) -= 2.0*coll_deriv.col(i).transpose()*v;
         //grad_eigen(i) += (invcov*dMtot[i]).trace();
      }

      grad[0]     = grad_eigen(0);
      grad[1]     = grad_eigen(1);
      grad[2]     = grad_eigen(2);

      //std::cout << "fun " << fun << std::endl;
      //std::cout << "dfx, dfy, dfz  = " << dfx << ", " << dfy << ", " << dfz << std::endl;
      //std::cout << "grad_eigen(0), grad_eigen(1), grad_eigen(2)  = " << grad_eigen(0) << ", " << grad_eigen(1) << ", " << grad_eigen(2) << std::endl;
      return fun;
    }

  };
}
*/
//new implementation using MFA 
inline FitResult coreFC( Block<real_t>* b, //replace with block that has encoded MFA model
			 diy::Master::ProxyWithLink const& cp,
			 Eigen::VectorXd const & fake_data,
			 Eigen::VectorXd const & v_coll, 
			 Eigen::VectorXd const & full, 
			 int dim1,
			 int dim2,
			 int dim3,
			 int i_grid,
			 int o_grid,
			 Eigen::MatrixXd const & mat_frac,
			 Eigen::MatrixXd const & mat_trans,
			 Eigen::MatrixXd const & mat_dirt,
			 Eigen::MatrixXd const & mat_data,
			 std::ofstream& debugFile,
			 bool doScan,
			 bool pearson,
			 size_t max_number_iterations,
			 double factr,
			 double pgtol
			)
{
  auto start = high_resolution_clock::now();
  verbose=false;
  bool write=false;
  float last_chi_min = FLT_MAX;
  int best_iter = -1;
  float global_chi_min = FLT_MAX;
  int best_grid_point = -99;
  size_t n_iter = 0;
  float chi_min = FLT_MAX;
  double current_best_grid_point;
  double best_grid_pointx;
  double best_grid_pointy;
  double best_grid_pointz;
  int flag = -99;
  int best_flag = -99;
  unsigned int num_iters = 0;
  unsigned int num_fun_calls = 0;
  
  Eigen::MatrixXd invcovcol;
  Eigen::MatrixXd covcol;

  double t1 = MPI_Wtime(); 
  if( !pearson ) {
    auto const & rescnp = updateInvCovCNP(mat_frac, mat_trans, mat_dirt, mat_data, full, fake_data, debugFile);
    invcovcol = std::get<1>(rescnp);
    covcol    = std::get<0>(rescnp);
  } else {
    auto const & rescnp = updateInvCov(mat_frac, mat_trans, mat_dirt, mat_data, full);
    invcovcol = std::get<1>(rescnp);
    covcol    = std::get<0>(rescnp);
  }
  double t2 = MPI_Wtime();
  debugFile << "fake_data: [" << fake_data.transpose() << "]" << std::endl;
  debugFile << "oscillated spectrum: [" << v_coll.transpose() << "]" << std::endl;
  debugFile << "fake_data - osc spectrum: [" << (fake_data - v_coll).transpose() << "]" << std::endl;
  debugFile << "inverse collapsed matrix: [" << invcovcol.diagonal().transpose() << "]" << std::endl;
  //std::cout << "diff: " << (fake_data - v_coll) << std::endl;
  double this_chi = calcChi(fake_data - v_coll, invcovcol);
  debugFile << "this_chi: " << this_chi << std::endl;
  double t3 = MPI_Wtime();
    
  //if(this_chi<0) 
  std::cout << "i_grid (dmsq, sin2theta_mue, sin2theta24), this_chi: " << i_grid << "(" << int(i_grid/dim2/dim3) << ", " << int(i_grid/dim2)%dim3 << ", " << i_grid%dim3 << "), " << this_chi << std::endl; 
 
  //do the optimization if doScan flag == 0
  bool debug = false;
  /*if( !doScan ){
    // //decode the grid point in 3d
    int gridx_index = int(i_grid/(dim2*dim3));
    int gridy_index = int(i_grid/dim2)%dim3;
    int gridz_index = int(i_grid%dim2);
    
    if(verbose) 
      std::cout << "i_grid, x, y, z = " << i_grid << ", " << gridx_index << ", " << gridy_index << ", " << gridz_index << std::endl;
    
    std::array<double, 3> x0{double(gridx_index), double(gridy_index), double(gridz_index)};
    const std::array<double, 3> lb{0., 0., 0.};
    const std::array<double, 3> ub{dim1-1.0, dim2-1.0, dim3-1.0};
    // 0 if unbounded,
    // 1 if only a lower bound,
    // 2 if both lower and upper bounds,
    // 3 if only an upper bound.
    const std::array<int, 3> bound_type{2, 2, 2};
    
    //pass the mfa model here which is already encode in the block
    auto startcputime = clock(); auto wcts = std::chrono::system_clock::now();
    if(verbose) std::cout << "max_number_iterations = " << max_number_iterations << std::endl;

    double step=1.;
    
    int n=3;
    std::mt19937 mt;
    mt.seed(0);
    std::random_device r;
    default_random_engine eng{r()};
    //debugFile << "n_iter, dim1, sinsq, chi_min, result.f_tol, result.num_iters, result.num_fun_calls, result.time_on_cauchy_point, result.time_on_line_search, result.time_on_subspace_minimization, flag \n";
    uniform_real_distribution<double> urd(0., 1.);
    while( n_iter<max_number_iterations){
      if(n_iter==0){
        x0[0]=gridx_index;
        x0[1]=gridy_index;
        x0[2]=gridz_index;
      } else {
        x0[0]=int(urd(mt)*(dim1));
        x0[1]=int(urd(mt)*dim2);
        x0[2]=int(urd(mt)*dim3);
      }
      n_iter++;
      //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
      lbfgsb::LLR fun(*b, cp, fake_data, mat_frac, mat_trans, mat_dirt, mat_data, dim1, dim2, dim3, pearson, n);
      
      lbfgsb::Optimizer optimizer{lb.size()};
      // Can adjust many optimization configs.
      // E.g. `iprint`, `factr`, `pgtol`, `max_iter`, `max_fun`, `time_limit_sec`
      //optimizer.iprint = -1;
      optimizer.factr = factr;
      optimizer.pgtol = pgtol;

      std::array<double, 3> grad;
      auto result = optimizer.minimize(fun, x0, lb.data(), ub.data(), bound_type.data());
      flag = result.warn_flag;
      if(result.task.find("NORM_OF_PROJECTED_GRADIENT") != std::string::npos ) flag = 4;
      if(result.task.find("REL_REDUCTION_OF_F_") != std::string::npos ) flag = 5;
      double chi_min = result.f_opt;
      if(chi_min < last_chi_min){
        //std::cout << "flag = " << flag << std::endl;
        last_chi_min = chi_min;
        best_grid_point  = current_best_grid_point;
        best_iter = n_iter;
	best_flag = flag;
        best_grid_pointx  = x0[0];
        best_grid_pointy  = x0[1];
        best_grid_pointz  = x0[2];
      }
      //debugFile << n_iter << ", " << roundoff(x0[0],2) << ", " << roundoff(x0[1],2) << ", " << roundoff(result.f_opt,3) << ", " << result.f_tol << ", " << result.num_iters << ", " << result.num_fun_calls << ", " << result.time_on_cauchy_points << ", " << result.time_on_line_search << ", " << result.time_on_subspace_minimization << ", " << flag << "\n";
    } // End loop over iterations

    global_chi_min = last_chi_min; //the global minimum chi2
    
    auto endcputime = clock(); auto wcte = std::chrono::system_clock::now();
    std::chrono::duration<double> wctduration = (wcte - wcts);
    int best_grid_point = x0[0]*dim2*dim3 + x0[1]*dim3 + x0[2];

  }// do minimization
  else*/ 
  last_chi_min = 0.;

  //assert that fakedataC should have the same dimension as collspec
  if(fake_data.size() != v_coll.size() ) std::cout << "check the collapsing method!" << std::endl;
  std::vector<double> fakedataC, collspec;
  for(uint i=0; i < fake_data.size(); i++) fakedataC.push_back(fake_data(i));
  for(uint i=0; i < v_coll.size(); i++) collspec.push_back(v_coll(i));
  
  FitResult fr = {best_iter, best_grid_point, best_grid_pointx, best_grid_pointy, best_grid_pointz, last_chi_min, this_chi-last_chi_min, best_flag};
  
  return fr;
}



void doScan(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
	    Eigen::MatrixXf const & mat_frac_, Eigen::MatrixXf const & mat_trans_,
	    Eigen::MatrixXf const & mat_dirt_, Eigen::MatrixXf const & mat_data_,
	    Eigen::VectorXd const & corecoll,  int dim1, int dim2, int dim3,
	    HighFive::File* file, std::vector<size_t> const & rankwork, 
	    double tol, double factr, size_t iter, bool pearson,
	    std::vector<float> const & xcoord, std::vector<float> const & ycoord, std::vector<float> const & zcoord, 
	    std::ofstream& debugFile, bool debug)
{

  // Cast MatrixXf to MatrixXd
  Eigen::MatrixXd mat_frac = mat_frac_.cast<double>();
  Eigen::MatrixXd mat_trans = mat_trans_.cast<double>();
  Eigen::MatrixXd mat_dirt = mat_dirt_.cast<double>();
  Eigen::MatrixXd mat_data = mat_data_.cast<double>();


  int nUniverses=1;
  double starttime, endtime;
  std::vector<FitResult> results;
  std::vector<int> v_grid, v_univ, v_iter, v_best;
  std::vector<double> v_last, v_dchi;
  std::vector<double> v_gridx, v_gridy, v_gridz;
  std::vector<std::vector<double> > v_fakedataC, v_collspec;
  std::vector<double> v_timeIteration, v_timeOptimizer, v_timeDecode, v_timeGrad;
  
  size_t pStart = rankwork[0];
  size_t uStart = rankwork[1];
  size_t pLast  = rankwork[2];
  size_t uLast  = rankwork[3];

  size_t i_begin = pStart * nUniverses + uStart;
  size_t i_end   = pLast  * nUniverses + uLast;

  //f/mt::print(stderr, "[{}] a,b,c,d: {} {} {} {} start at {} end at {}  lends {}\n", rank, pStart, uStart, pLast, uLast, i_begin, i_end, i_end-i_begin);
  size_t lenDS = i_end - i_begin;
    
  results.reserve(lenDS);
  v_grid.reserve(lenDS);
  v_univ.reserve(lenDS);
  v_iter.reserve(lenDS);
  v_best.reserve(lenDS);
  v_last.reserve(lenDS);
  v_dchi.reserve(lenDS);
  v_gridx.reserve(lenDS);
  v_gridy.reserve(lenDS);
  v_gridz.reserve(lenDS);
  v_fakedataC.reserve(lenDS);
  v_collspec.reserve(lenDS);
 
  system_clock::time_point t_init = system_clock::now();
  
    std::vector<std::vector<size_t> > RW;
    if (pStart == pLast) RW.push_back({pStart, uStart, uLast});
    else {
      RW.push_back({pStart, uStart, nUniverses});
       for (size_t _p = pStart+1; _p<pLast;++_p) {
          RW.push_back({_p, 0, nUniverses});
       }
       if (uLast>0) RW.push_back({pLast, 0, uLast});

    }

  for (auto r : RW) {
    
    int i_grid = r[0];
    if (debug && i_grid!=0) return;
    int x = int(i_grid/dim2/dim3);
    int y = int(i_grid/dim2)%dim3;
    int z = i_grid%dim3;
    //if(i_grid > 1000) continue;
    //if( z < y ) continue;
    //if( i_grid != 0 && i_grid != 30 && i_grid != 31 && ( int(i_grid/dim2/dim3) != 54 || int(i_grid/dim2)%dim3 != 3 || i_grid%dim3 != 31 ) ) continue;
    auto const & specfull = getSpectrum(b, cp, i_grid, dim2, dim3);
    //for(int bin = 0; bin < 1; bin++ ) std::cout << "i_grid, bin, spectrum: " << i_grid << ", " << bin << ", " << std::setprecision(8) << float(specfull(bin)) << std::endl; 
    double t1 = MPI_Wtime();
    auto const & speccoll = collapseVectorTrans(specfull, mat_trans);
    double t2 = MPI_Wtime();
    if(rank == 0) std::cout << "time to collapse spectrum: " << t2-t1 << " seconds" << std::endl;
    //auto const & corefull = getSpectrum(b, cp, 0, dim2, dim3); //test using null spectrum
    //auto const & corecoll = collapseVectorTrans(corefull, mat_trans);
    v_gridx.push_back(xcoord[x]);
    v_gridy.push_back(ycoord[y]);
    v_gridz.push_back(zcoord[z]);
    //std::cout <<  "*********** grid index, delta-msq, sin^2theta_{mue}, sin^2theta_{24} : " << i_grid << ", " << xcoord[x] << ", " << ycoord[y] << ", " << zcoord[z] << std::endl;
    debugFile << "\n *********** grid index, delta-msq, sin^2theta_{mue}, sin^2theta_{24} : " << i_grid << ", " << xcoord[x] << ", " << ycoord[y] << ", " << zcoord[z] << "\n" << std::endl;
    //debugFile << "corefull: [" << corefull.transpose() << "]" << std::endl;
    debugFile << "corecoll: [" << corecoll.transpose() << "]" << std::endl;
    debugFile << "specfull: [" << specfull.transpose() << "]" << std::endl;
    debugFile << "speccoll: [" << speccoll.transpose() << "]" << std::endl;
    
    starttime = MPI_Wtime();
    
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    //if(rank==0) std::cout << "rank, Memory usage beore coreFC: " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    double t3 = MPI_Wtime();
    results.push_back(coreFC(b, cp, corecoll, speccoll, specfull, dim1, dim2, dim3, i_grid, -1, mat_frac, mat_trans, mat_dirt, mat_data, debugFile, true, pearson, iter, factr, tol )); 
    double t4 = MPI_Wtime();
    getrusage(RUSAGE_SELF, &usage);
    //if(rank==0) std::cout << "rank, Memory usage after coreFC: " << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    //if(rank==0) std::cout << "time to coreFC: " << t4-t3 << " seconds" << std::endl;
    
    v_univ.push_back(0);
    v_grid.push_back(i_grid);
    
    endtime   = MPI_Wtime();
    system_clock::time_point now = system_clock::now();
    
    auto t_elapsed = now - t_init;
    auto t_togo = t_elapsed * (int(lenDS) - i_grid)/(i_grid+1);
    auto t_eta = now + t_togo;
    std::time_t t = system_clock::to_time_t(t_eta);
    
    if (rank==0 && i_grid%100==0) fmt::print(stderr, "[{}] gridp {}/{} took {} seconds. ETA: {}",cp.gid(), i_grid, lenDS, endtime-starttime, std::ctime(&t));
  }
  
  // Write to HDF5
  starttime   = MPI_Wtime();
  HighFive::DataSet d_last_chi_min    = file->getDataSet("last_chi_min"   );
  HighFive::DataSet d_delta_chi       = file->getDataSet("delta_chi"      );
  HighFive::DataSet d_best_grid_point = file->getDataSet("best_grid_point");
  HighFive::DataSet d_n_iter          = file->getDataSet("n_iter"         );
  // write out this grid and universe
  HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
  HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
  HighFive::DataSet d_gridx          = file->getDataSet("gridx");
  HighFive::DataSet d_gridy          = file->getDataSet("gridy");
  HighFive::DataSet d_gridz          = file->getDataSet("gridz");
  
  size_t d_bgn = rankwork[0];
  for (auto res : results) {
    v_iter.push_back(res.n_iter);
    v_best.push_back(res.best_grid_point);
    v_last.push_back(res.last_chi_min);
    v_dchi.push_back(res.delta_chi);
  }
 

  double minchi = *min_element(v_dchi.begin(), v_dchi.end());
  std::map<double,int> mapbp;
  //std::cout << "i_grid, minchi: " << i_grid << ", " << minchi << std::endl;
  for (int i=0; i < v_dchi.size(); i++ ){
      v_last[i] = minchi;
      double this_chi = v_dchi[i];
      mapbp[this_chi] = i;
      v_dchi[i] -= minchi;
      if(v_dchi[i] == 0) std::cout << "grid, this_chi, minchi, deltachi: " << i << ", " << this_chi << ", " << v_last[i] << ", " << v_dchi[i] << std::endl;
      //std::cout << i << ", " << this_chi << ", " << v_last[i] << ", " << v_dchi[i] << std::endl;
  }

  //find best point
  auto it = mapbp.find(minchi);

  if (it != mapbp.end()) {
      std::cout << "The integer associated with " << minchi << " is " << it->second << std::endl;
  }

  d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
  d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
  d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
  d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
  d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
  d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
  d_gridx.select(            {d_bgn, 0}, {size_t(v_gridx.size()), 1}).write(v_gridx);
  d_gridy.select(            {d_bgn, 0}, {size_t(v_gridy.size()), 1}).write(v_gridy);
  d_gridz.select(            {d_bgn, 0}, {size_t(v_gridz.size()), 1}).write(v_gridz);
  if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
}

/*
void doMin(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
	   const char * xmldata, sbn::SBNconfig const & myconf,
	   TMatrixD const & covmat, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
	   SignalGenerator signal,
	   HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
	   double tol, size_t iter, bool debug, bool noWrite=false, int msg_every=100)
{ 
  
  //return 0;
}

// TODO add size_t writeEvery to prevent memory overload
void doFC(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
	  const char * xmldata, sbn::SBNconfig const & myconf, int mass, int dim2, int dim3,
	  Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG, Eigen::VectorXd const & ecore,
	  SignalGenerator signal, HighFive::File* file, ofstream &chi2value, ofstream &debugFile,  
	  std::vector<size_t> const & rankwork, int nUniverses, int outer_gpt, int inner_gpt,
	  double tol, double factr, size_t iter, int degree, bool debug, bool doFullGrid, bool noWrite=false, int msg_every=1 )
{
  //std::cout << "doFC, signal.gridsize() " <<  signal.gridsize() << std::endl;
  double starttime, endtime;
  starttime = MPI_Wtime();
  double timeLMAT = 0.;
  double timeFC = 0.;
  double timeIteration = 0.;
  double timeOptimizer = 0.;
  double timeDecode = 0.;
  double timeGrad = 0.;
  std::vector<FitResult> results;
  std::vector<int> v_grid, v_univ, v_iter, v_best, v_flag;
  std::vector<unsigned int> v_num_iters;
  std::vector<unsigned int> v_num_fun_calls;
  std::vector<double> v_previous_fval;
  std::vector<double> v_f_tol;
  std::vector<double> v_bestx, v_besty;
  std::vector<double> v_time_on_cauchy_points;
  std::vector<double> v_time_on_subspace_minimization;
  std::vector<double> v_time_on_line_search;
  std::vector<int> v_max_iter_exceeded;
  std::vector<int> v_max_fun_exceeded;
  std::vector<int> v_time_limit_exceeded;
  std::vector<double> v_last, v_dchi;
  std::vector<double> v_timeLMAT, v_timeFC, v_timeFCfd, v_timeIteration, v_timeOptimizer, v_timeDecode, v_timeGrad;
  std::vector<std::vector<double> > v_fakedataC, v_collspec;//, v_outpt, v_chi2vec, v_specbestC, v_invcovbestC;
  size_t pStart = rankwork[0];
  size_t uStart = rankwork[1];
  size_t pLast  = rankwork[2];
  size_t uLast  = rankwork[3];

  size_t i_begin = pStart * nUniverses + uStart;
  size_t i_end   = pLast  * nUniverses + uLast;

  //fmt::print(stderr, "[{}] a,b,c,d: {} {} {} {} start at {} end at {}  lends {}\n", rank, pStart, uStart, pLast, uLast, i_begin, i_end, i_end-i_begin);
  size_t lenDS = i_end - i_begin;

  if (!noWrite) {
    results.reserve(lenDS);
    v_grid.reserve(lenDS);
    v_univ.reserve(lenDS);
    v_iter.reserve(lenDS);
    v_best.reserve(lenDS);
    v_bestx.reserve(lenDS);
    v_besty.reserve(lenDS);
    v_last.reserve(lenDS);
    v_dchi.reserve(lenDS);
    v_flag.reserve(lenDS);
    v_num_iters.reserve(lenDS);
    v_num_fun_calls.reserve(lenDS);
    v_previous_fval.reserve(lenDS);
    v_f_tol.reserve(lenDS);
    v_time_on_cauchy_points.reserve(lenDS);
    v_time_on_subspace_minimization.reserve(lenDS);
    v_time_on_line_search.reserve(lenDS);
    v_max_iter_exceeded.reserve(lenDS);
    v_max_fun_exceeded.reserve(lenDS);
    v_time_limit_exceeded.reserve(lenDS);
    v_timeLMAT.reserve(lenDS);
    v_timeFC.reserve(lenDS);
    v_timeFCfd.reserve(lenDS);
    v_timeIteration.reserve(lenDS);
    v_timeOptimizer.reserve(lenDS);
    v_timeDecode.reserve(lenDS);
    v_timeGrad.reserve(lenDS);
    v_fakedataC.reserve(lenDS);
    v_collspec.reserve(lenDS);
  }

  //std::cout << "doFC" << std::endl;

  std::vector<std::vector<size_t> > RW;
  if (pStart == pLast) RW.push_back({pStart, uStart, uLast});
  else {
    RW.push_back({pStart, uStart, nUniverses});
     for (size_t _p = pStart+1; _p<pLast;++_p) {
        RW.push_back({_p, 0, nUniverses});
     }
     if (uLast>0) RW.push_back({pLast, 0, uLast});

  }

  system_clock::time_point t_init = system_clock::now();
  auto startcputime = clock(); auto wcts = std::chrono::system_clock::now();
  // write array of outer grid here (always 25x25x25)
  //float masmin = -1.995;
  //float ue4min = -2.0;
  //float umu4min = -2.0;
  //std::vector<float> dmasssqs, ue4s, umu4s;

  //for (size_t w=0; w < 25; w++) dmasssqs.push_back( roundoff(masmin+w*0.16,3));
  //for (size_t w=0; w < 25; w++) ue4s.push_back( roundoff(ue4min+w*0.06795937566,2));// std::cout << "ue4s.back(): " << ue4s.back() << std::endl;}
  //for (size_t w=0; w < 25; w++) umu4s.push_back( roundoff(umu4min+w*0.06795937566,2));
  
  for (auto r : RW) {

     size_t i_grid = r[0];
     if ( !doFullGrid && i_grid != outer_gpt ) continue;
     inner_gpt = transformtoinnergpt(i_grid);
     outer_gpt = transformtooutergpt(inner_gpt);
     //std::cout << "2. i_grid, inner_gpt, outer_gpt: " << i_grid << ", " << inner_gpt << ", " << outer_gpt << std::endl;

     //GridPoint gp = signal.getgrid(i_grid);
     //if(std::find(dmasssqs.begin(), dmasssqs.end(), roundoff(log10(gp[0]*gp[0]),3)) == dmasssqs.end()) continue;
     //if(std::find(ue4s.begin(), ue4s.end(), roundoff(log10(gp[1]),2)) == ue4s.end()) continue;
     //if(std::find(umu4s.begin(), umu4s.end(), roundoff(log10(gp[2]),2)) == umu4s.end()) continue;
 
     //std::cout << "outer_gpt = " << outer_gpt << std::endl;
     auto const & specfull = signal.predict(inner_gpt, false);
     auto const & speccoll   = collapseVectorEigen(specfull, myconf);
     auto const & corecoll   = collapseVectorEigen(ecore, myconf);

     if(verbose) 
       std::cout << "corecoll: " << corecoll << std::endl;
     //debugFile << "corecoll: " << corecoll << "\n";
     if(verbose) 
       std::cout << "specfull: " << specfull << std::endl;
     //debugFile << "specfull: " << specfull << "\n";
     if(verbose) 
       std::cout << "speccoll: " << speccoll << std::endl;
     //debugFile << "speccoll: " << speccoll << "\n";

     double starttimeLMAT = MPI_Wtime();
     Eigen::MatrixXd const & LMAT = cholDcollapsed(ECOV, specfull, myconf);
     double endtimeLMAT = MPI_Wtime();
     
     //creates a linear chain of blocks like in https://urldefense.proofpoint.com/v2/url?u=https-3A__github.com_diatomic_diy_blob_master_examples_simple_iexchange-2Dparticles.cpp-3F&d=DwIGAg&c=gRgGjJ3BkIsb5y6s49QqsA&r=k6h-opkRyS6dcFGVdbiNygHT00lxVm2iqC2QgLizjsQ&m=ADWd_Es8iKxWTqNtSAGTb63W8UbsqxA-jyGqZJ-fvH8&s=kzEBHQze8RX90NqbZnC5wxn1mbq1TmpHeK4KMabOqco&e= 
     for (int uu=r[1]; uu<r[2];++uu) {
       std::mt19937 rng((i_grid+1)*(uu+1)); // Mersenne twister

       double starttimeFC = MPI_Wtime();
       //Katie sampled directly from collapsed matrix
       //auto const & fake_data = poisson_fluctuate(sample(specfull, LMAT, rng),rng);//
       auto const & fake_dataC = poisson_fluctuate(sample(speccoll, LMAT, rng),rng);
       if( verbose ) 
         std::cout << "fake_dataC: " << fake_dataC.transpose() << std::endl;
       //debugFile << "fake_dataC: " << fake_dataC << "\n";
       double endtimeFCfd = MPI_Wtime();

       //std::cout << "coreFC" << std::endl;
       results.push_back(coreFC(fake_dataC, speccoll, specfull, b, cp, signal, mass, dim2, dim3, inner_gpt, i_grid, INVCOVBG, ECOV, timeIteration, timeOptimizer, timeDecode, timeGrad, myconf, debugFile, false, iter, factr, tol)); 
       double endtimeFC = MPI_Wtime();
       v_univ.push_back(uu);
       v_grid.push_back(outer_gpt);
       v_timeFC.push_back(endtimeFC-starttimeFC);
       v_timeFCfd.push_back(endtimeFCfd-starttimeFC);
       v_timeIteration.push_back(timeIteration);
       v_timeOptimizer.push_back(timeOptimizer);
       v_timeDecode.push_back(timeDecode);
       v_timeGrad.push_back(timeGrad);
       v_timeLMAT.push_back(endtimeLMAT-starttimeLMAT);
     }
     endtime   = MPI_Wtime();
     system_clock::time_point now = system_clock::now();
     
     auto endcputime = clock(); auto wcte = std::chrono::system_clock::now();
     std::chrono::duration<double> wctduration = (wcte - wcts);

     auto t_elapsed = now - t_init;
     auto t_togo = t_elapsed * (int(lenDS) - results.size())/(results.size()+1);
     auto t_eta = now + t_togo;
     std::time_t t = system_clock::to_time_t(t_eta);
     
     if (cp.gid()==0 && results.size()%msg_every==0) fmt::print(stderr, "[{}] gridp {}/{} ({} universes) took {} seconds. ETA: {}",cp.gid(), results.size(), lenDS, nUniverses, endtime-starttime, std::ctime(&t));
  }
  
  if (!noWrite) {
    
    // Write to HDF5
    starttime   = MPI_Wtime();
    HighFive::DataSet d_last_chi_min    = file->getDataSet("last_chi_min"   );
    HighFive::DataSet d_delta_chi       = file->getDataSet("delta_chi"      );
    HighFive::DataSet d_best_grid_point = file->getDataSet("best_grid_point");
    HighFive::DataSet d_best_grid_point_x = file->getDataSet("best_grid_point_x");
    HighFive::DataSet d_best_grid_point_y = file->getDataSet("best_grid_point_y");
    HighFive::DataSet d_n_iter          = file->getDataSet("n_iter"         );
    HighFive::DataSet d_warn_flag       = file->getDataSet("warn_flag"    );
    HighFive::DataSet d_num_iters       = file->getDataSet("num_iters");
    HighFive::DataSet d_num_fun_calls   = file->getDataSet("num_fun_calls");
    HighFive::DataSet d_previous_fval   = file->getDataSet("previous_fval");
    HighFive::DataSet d_f_tol           = file->getDataSet("f_tol");
    HighFive::DataSet d_time_on_cauchy_points = file->getDataSet("time_on_cauchy_points");
    HighFive::DataSet d_time_on_subspace_minimization = file->getDataSet("time_on_subspace_minimization");
    HighFive::DataSet d_time_on_line_search = file->getDataSet("time_on_line_search");
    HighFive::DataSet d_max_iter_exceeded = file->getDataSet("max_iter_exceeded");
    HighFive::DataSet d_max_fun_exceeded = file->getDataSet("max_fun_exceeded");
    HighFive::DataSet d_time_limit_exceeded = file->getDataSet("time_limit_exceeded");
    // write out this grid and universe
    HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
    HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
    HighFive::DataSet d_timeLMAT          = file->getDataSet("timeLMAT");
    HighFive::DataSet d_timeFC          = file->getDataSet("timeFC");
    HighFive::DataSet d_timeFCfd          = file->getDataSet("timeFCfd");
    HighFive::DataSet d_timeIteration          = file->getDataSet("timeIteration");
    HighFive::DataSet d_timeOptimizer          = file->getDataSet("timeOptimizer");
    HighFive::DataSet d_timeDecode          = file->getDataSet("timeDecode");
    HighFive::DataSet d_timeGrad          = file->getDataSet("timeGrad");
    
    size_t d_bgn = i_begin;//rankwork[0]*nUniverses;
    for (auto res : results) {
      v_iter.push_back(res.n_iter);
      v_best.push_back(res.best_grid_point);
      v_bestx.push_back(res.best_grid_point_x);
      v_besty.push_back(res.best_grid_point_y);
      v_last.push_back(res.last_chi_min);
      v_flag.push_back(res.flag);
      v_num_iters.push_back(res.num_iters);
      v_num_fun_calls.push_back(res.num_fun_calls);
      v_previous_fval.push_back(res.previous_fval);
      v_f_tol.push_back(res.f_tol);
      v_time_on_cauchy_points.push_back(res.time_on_cauchy_points);
      v_time_on_subspace_minimization.push_back(res.time_on_subspace_minimization);
      v_time_on_line_search.push_back(res.time_on_line_search);
      v_max_iter_exceeded.push_back(int(res.max_iter_exceeded));
      v_max_fun_exceeded.push_back(int(res.max_fun_exceeded));
      v_time_limit_exceeded.push_back(int(res.time_limit_exceeded));
      v_dchi.push_back(res.delta_chi);
    }
    int lenVec=v_grid.size()-1;
    d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
    d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
    d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
    d_best_grid_point_x.select(  {d_bgn, 0}, {size_t(v_bestx.size()), 1}).write(v_bestx);
    d_best_grid_point_y.select(  {d_bgn, 0}, {size_t(v_besty.size()), 1}).write(v_besty);
    d_warn_flag.select(  {d_bgn, 0}, {size_t(v_flag.size()), 1}).write(v_flag);
    d_num_iters.select(  {d_bgn, 0}, {size_t(v_num_iters.size()), 1}).write(v_num_iters);
    d_num_fun_calls.select(  {d_bgn, 0}, {size_t(v_num_fun_calls.size()), 1}).write(v_num_fun_calls);
    d_previous_fval.select(  {d_bgn, 0}, {size_t(v_previous_fval.size()), 1}).write(v_previous_fval);
    d_f_tol.select(  {d_bgn, 0}, {size_t(v_f_tol.size()), 1}).write(v_f_tol);
    d_time_on_cauchy_points.select(  {d_bgn, 0}, {size_t(v_time_on_cauchy_points.size()), 1}).write(v_time_on_cauchy_points);
    d_time_on_subspace_minimization.select(  {d_bgn, 0}, {size_t(v_time_on_subspace_minimization.size()), 1}).write(v_time_on_subspace_minimization);
    d_time_on_line_search.select(  {d_bgn, 0}, {size_t(v_time_on_line_search.size()), 1}).write(v_time_on_line_search);
    d_max_iter_exceeded.select(  {d_bgn, 0}, {size_t(v_max_iter_exceeded.size()), 1}).write(v_max_iter_exceeded);
    d_max_fun_exceeded.select(  {d_bgn, 0}, {size_t(v_max_fun_exceeded.size()), 1}).write(v_max_fun_exceeded);
    d_time_limit_exceeded.select(  {d_bgn, 0}, {size_t(v_time_limit_exceeded.size()), 1}).write(v_time_limit_exceeded);
    d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
    d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
    d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
    d_timeLMAT.select(           {d_bgn, 0}, {size_t(v_timeLMAT.size()), 1}).write(v_timeLMAT);
    d_timeFC.select(           {d_bgn, 0}, {size_t(v_timeFC.size()), 1}).write(v_timeFC);
    d_timeFCfd.select(           {d_bgn, 0}, {size_t(v_timeFCfd.size()), 1}).write(v_timeFCfd);
    d_timeIteration.select(      {d_bgn, 0}, {size_t(v_timeIteration.size()), 1}).write(v_timeIteration);
    d_timeOptimizer.select(      {d_bgn, 0}, {size_t(v_timeOptimizer.size()), 1}).write(v_timeOptimizer);
    d_timeDecode.select(         {d_bgn, 0}, {size_t(v_timeDecode.size()), 1}).write(v_timeDecode);
    d_timeGrad.select(           {d_bgn, 0}, {size_t(v_timeGrad.size()), 1}).write(v_timeGrad);
    endtime   = MPI_Wtime();
    if( cp.gid() == 0 ) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
  }
}
*/
bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// --- main program ---//
int main(int argc, char* argv[]) {
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    double T0   = MPI_Wtime();

    time_t now;
    time (&now);
    if (world.rank()==0) fmt::print(stderr, "Start at {}", std::ctime(&now));

    size_t nPoints =-1;
    int mode=0;
    int msg_every=100;
    size_t nUniverses=1;
    int NTEST(0);
    int mfadim=3;
    std::string out_file="test.h5";
    //std::string time_file="timestamps.h5";
    std::string tag="";
    //std::string xml="";
    std::string infile  = "";                 // h5 input file
    std::string fileName  = "";                 // root input file
    std::string debugFileName = "";                 // debug file name

    double tol(1.000e-4);
    double factr(1.000e+7);
    size_t iter(20);
    int degree = 2;
    int fgpt(0);
    int lgpt(0);

    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    //ops >> Option("timefile",     time_file,   "Output filename for timestamps.");
    ops >> Option('i', "input",      infile,     "Input filename.");
    ops >> Option('o', "output",     out_file,   "Output filename.");
    ops >> Option('u', "nuniverse",  nUniverses, "Number of universes");
    ops >> Option("ntest",           NTEST ,     "Test Run");
    ops >> Option("tol",             tol,        "Minimiser gradient tolerance");
    ops >> Option("factr",           factr,      "Minimiser factr");
    ops >> Option("iter",            iter,       "Max number of iterations.");
    ops >> Option('t', "tag",        tag,        "Tag.");
    ops >> Option('f', "rootfile",   fileName,   "root file");
    ops >> Option('l', "log file",   debugFileName,   "debugging logfile");
    ops >> Option("degree",          degree,     "MFA science degree");
    ops >> Option("fgpt",	     fgpt,       "First grid point");
    ops >> Option("lgpt",	     lgpt,       "Last grid point");
    ops >> Option("mfadim",          mfadim,     "Dimension of the MFA model");
    ops >> Option("msg",             msg_every,  "Print a progress message every m gridpoints on rank 0 to stderr.");
    ops >> Option("mode",            mode, "Mode 0 is default --- dimension 2 is electron, mode 1 is muon");
    bool debug      = ops >> Present('d', "debug", "Operate on single gridpoint only");
    bool doFullGrid = ops >> Present("full", "Do full grid instead of single outer gridpoint");
    bool statonly   = ops >> Present("stat", "Statistical errors only");
    bool nowrite    = ops >> Present("nowrite", "Don't write output --- for performance estimates only");
    bool simplescan = ops >> Present("scan", "Simple scan, no FC");
    bool pearson    = ops >> Present("pearson", "Use Pearson chi2");
    if (ops >> Present('h', "help", "Show help")) {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }
    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running SBN Feldman Cousins ***\n");
    }

    // Whole bunch of tests
    if ( world.rank() == 0 ) {
       //std::cout << "world.size() = " << world.size() << std::endl;
       if (int(world.size()) > nPoints) {
          std::cerr << "Impossible to run on more ranks than grid points, exiting.\n";
          exit(1);
       }
       std::vector<std::string> infiles = {};
       for (auto f : infiles) {
          if (!file_exists(f)) { 
             std::cerr << "Specified input file " << f <<" does not exist, exiting\n";
             exit(1);
          }
       }
       if (tag=="") {
          std::cerr << "tag (-t, --tag) cannot be undefined, exiting\n";
          exit(1);
       }
    }


    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage 0a: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    double t1 = MPI_Wtime();
    std::vector<float> v_buff;
    std::vector<double> v_buffd;
    int ncols, nrows;
 

    // Open the ROOT file
    std::string matrixName = "matrix_data_total";
     
    TFile* file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open the file '" << fileName << "'" << std::endl;
        return 1;
    }

    // Get the TMatrixT object from the file
    TMatrixT<Double_t>* matrix = dynamic_cast<TMatrixT<Double_t>*>(file->Get(matrixName.c_str()));
    if (!matrix) {
        std::cerr << "Failed to retrieve the matrix '" << matrixName << "' from the file" << std::endl;
        file->Close();
        return 1;
    }

    // Get the number of bins in the z axis
    Int_t numBinsZ = matrix->GetNrows();
    Int_t numBinsX = matrix->GetNcols();
    std::vector<double> _core;
    // Extract the specified column as a TMatrixT object
    TMatrixT<Double_t> columnMatrix = matrix->GetSub(0, numBinsZ - 1, 0, numBinsX - 1);
    // Iterate over the elements in the column
    for (Int_t i = 0; i < numBinsZ; i++) {
        for (Int_t j = 0; j < numBinsX; j++) {
            _core.push_back(columnMatrix(i, j));  // Access the element in the column
        }
    }
    diy::mpi::broadcast(world, _core, 0);

    Eigen::Map<Eigen::VectorXd> ecore(_core.data(), _core.size(), 1);
   
    if (world.rank() == 0) loadData(Form("%s",infile.c_str()), "matrix_fractional_flux_Xs_G4_det", v_buff,   nrows, ncols);
    Eigen::MatrixXf _mat_frac = bcMatrixXf(world, v_buff, nrows, ncols);
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage before matrix_fractional_flux_Xs_G4_det: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;

    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage after matrix_fractional_flux_Xs_G4_det: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
  
    if (world.rank() == 0) loadData(Form("%s",infile.c_str()), "matrix_transition", v_buff,   nrows, ncols);
    Eigen::MatrixXf _mat_trans = bcMatrixXf(world, v_buff, nrows, ncols);
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage matrix_transition: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
  
    if (world.rank() == 0) loadData(Form("%s",infile.c_str()), "matrix_data_1xN", v_buff,   nrows, ncols);
    Eigen::MatrixXf _mat_data = bcMatrixXf(world, v_buff, nrows, ncols);
    Eigen::VectorXd dataspec = _mat_data.row(0).cast<double>();
    std::cout << "dataspec size: " << dataspec.size() << std::endl;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage matrix_data_1xN: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
  
    if (world.rank() == 0) loadData(Form("%s",infile.c_str()), "matrix_absolute_dirt_MCstat", v_buff,   nrows, ncols);
    Eigen::MatrixXf _mat_dirt = bcMatrixXf(world, v_buff, nrows, ncols);
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rank, Memory usage matrix_absolute_dirt_MCstat: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
 
    double t2 = MPI_Wtime();
    //load spectrum at rank 0 
    std::vector<float> xcoord;
    std::vector<float> ycoord;
    std::vector<float> zcoord;
    int nBins = 0;
    //Eigen::Tensor<float,4> my_grid_;
    std::vector<float> myVector;
    //std::fill(my_grid.begin(), my_grid.end(), 0);

    getrusage(RUSAGE_SELF, &usage);
    if(world.rank()==0) std::cout << "rank, Memory usage before reading grid: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;

    double t3 = MPI_Wtime();
    if (world.rank()==0){
    	myVector = makeSignalObject(infile, xcoord, ycoord, zcoord, nBins);
    }
    double t4 = MPI_Wtime();
    if(world.rank() == 0) std::cout << "time to read input: " << t2-t1 << " seconds" << std::endl;
    if(world.rank() == 0) std::cout << "time to create tensor: " << t4-t3 << " seconds" << std::endl;
    getrusage(RUSAGE_SELF, &usage);
    if(world.rank()==0) std::cout << "rank, Memory usage after reading grid: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    diy::mpi::broadcast(world, xcoord, 0);   
    diy::mpi::broadcast(world, ycoord, 0);   
    diy::mpi::broadcast(world, zcoord, 0);   
    diy::mpi::broadcast(world, myVector, 0);  
    diy::mpi::broadcast(world, nBins, 0);   

    world.barrier();
    int dim1  = xcoord.size();
    int dim2  = ycoord.size();
    int dim3  = zcoord.size();
    if(world.rank()==0)  std::cout << "rank, dim1, dim2, dim3, nbins, myVector, size: " << world.rank() << ", " << dim1 << ", " << dim2 << ", " << dim3 << ", " << nBins << ", " << myVector.size() << std::endl;
    if(fgpt == 0 && lgpt == 0 ) nPoints = dim1*dim2*dim3;
    else nPoints = fabs(lgpt-fgpt);
    //initiate vector here:
    //std::vector<float> my_grid(nPoints*nBins);
    //transform vector back to Eigen::Tensor
    //declare the new eigen tensor
    getrusage(RUSAGE_SELF, &usage);
    if(world.rank()==0) std::cout << "rank, Memory usage before initilizing tensor: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;

    double time0 = MPI_Wtime();

    if( world.rank()==0 ) {
      fmt::print(stderr, "***********************************\n");
      fmt::print(stderr, "    Output will be written to {}\n", out_file);
      fmt::print(stderr, "    Inner Points:    {}\n"         , nPoints);
      fmt::print(stderr, "    Outer Points:    {}\n"         , nPoints);
      fmt::print(stderr, "    nBins:     {}\n"               , nBins);
      fmt::print(stderr, "    Universes: {}\n"               , nUniverses);
      fmt::print(stderr, "    Total size of dataset:  {}\n"  , nPoints*nUniverses);
      fmt::print(stderr, "    first gridpoint :  {}\n"       , fgpt);
      fmt::print(stderr, "    last gridpoint :  {}\n"        , lgpt);
      fmt::print(stderr, "    iter :  {}\n"                  , iter);
      fmt::print(stderr, "    tol:    {}\n"                  , tol);
      fmt::print(stderr, "    factr:  {}\n"                  , factr);
      if (statonly) fmt::print(stderr, "    S T A T  O N L Y \n");
      if (debug)    fmt::print(stderr,    "    D E B U G \n"    );
      if (nowrite)  fmt::print(stderr,  "  N O   W R I T E \n"  );
      if (simplescan)  fmt::print(stderr,  "       S C A N \n"  );
      if (pearson)  fmt::print(stderr, "     Pearson chi^2 \n"    );
      else          fmt::print(stderr, "         CNP chi^2 \n"    );
      fmt::print(stderr, "***********************************\n");
    }

    // Create hdf5 file structure here 
    HighFive::File* f_out  = new HighFive::File(out_file,
                        HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));

    // Create datasets needed TODO nbinsC --- can we get that from somewhere?
    createDataSets(f_out, nPoints, nUniverses);
   
    // First rank also writes the grid so we know what the poins actually are
    //if (world.rank() == 0)  writeGrid(f_out, myVector, xcoord, ycoord, zcoord);

    // Now more blocks as we have universes
    size_t blocks = world.size();//nPoints;// *nUniverses;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds<real_t> fc_domain(3);
    fc_domain.min[0] = 0.;
    fc_domain.max[0] = blocks-1;
    fc_domain.min[1] = 0.;
    fc_domain.max[1] = blocks-1;
    fc_domain.min[2] = 0.;
    fc_domain.max[2] = blocks-1;

    diy::FileStorage               storage("./DIY.XXXXXX");
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds<real_t>> fc_decomposer(mfadim, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, mfadim, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, mfadim, true);
    diy::Master                    fc_master(world, 1, -1, &Block<real_t>::create, &Block<real_t>::destroy, &storage, &Block<real_t>::save, &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    fc_decomposer.decompose(world.rank(),
                         assigner,
                         [&](int gid, const Bounds<real_t>& core, const Bounds<real_t>& bounds, const Bounds<real_t>& domain, const RCLink<real_t>& link)
                         { Block<real_t>::add(gid, core, bounds, domain, link, fc_master, mfadim, nBins+mfadim, 0.0); });

    double T10 = MPI_Wtime();
    getrusage(RUSAGE_SELF, &usage);
    if(world.rank()==0)  std::cout << "rank, Memory usage 4: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    fc_master.foreach([world, myVector, dim1, dim2, dim3, nBins, mfadim, degree](Block<real_t>* b, const diy::Master::ProxyWithLink& cp){
    	makeSignalModel(world, b, cp, myVector, dim1, dim2, dim3, nBins, mfadim, degree);});
    getrusage(RUSAGE_SELF, &usage);
    if(world.rank()==0) std::cout << "rank, Memory usage 5: " << world.rank() << ", " << usage.ru_maxrss/1000000. << " GB" << std::endl;
    double T11   = MPI_Wtime();

    if (world.rank()==0) std::cout << "time to build model: " << T11-T10 << " seconds." << std::endl;

    //do LoadBalance
    size_t _S(nPoints*nUniverses);
    std::vector<size_t> _L;
    size_t maxwork = size_t(ceil(_S/world.size()));
    //std::cout << "maxwork, world.size(): " << maxwork << ", " << world.size() << std::endl;
    for (size_t r=0; r <                _S%world.size(); ++r) _L.push_back(maxwork + 1);
    for (size_t r=0; r < world.size() - _S%world.size(); ++r) _L.push_back(maxwork    );


    std::vector<size_t> _bp, _bu;
    _bp.push_back(0);
    _bu.push_back(0);

    size_t _n(0), _temp(0);
    for (size_t i=0; i<nPoints;++i) {
       for (size_t j=0; j<nUniverses;++j) {
          if (_temp == _L[_n]) {
             _bp.push_back(i);
             _bu.push_back(j);
             _temp = 0;
             _n+=1;
          }
          _temp+=1;
       }
    }
    _bp.push_back(nPoints-1);
    _bu.push_back(nUniverses);
    
    std::vector<size_t> rankwork;
    fmt::print(stderr, "[{}] Got {} ranks and {} sets of work {}\n", world.rank(), world.size(), _bu.size() -1, _bp.size() -1);
    world.barrier();
    rankwork.push_back(_bp[world.rank()]);
    rankwork.push_back(_bu[world.rank()]);
    rankwork.push_back(_bp[world.rank()+1]);
    rankwork.push_back(_bu[world.rank()+1]);

    double T12 = MPI_Wtime();
    double starttime = MPI_Wtime();
  
    //declare the offstream object to write matrices and spectrums	
    std::ofstream debugFile; 
    debugFile.open(debugFileName.c_str());

    // Check if the file was successfully opened
    if (!debugFile.is_open()) {
        // Handle the error (e.g., print an error message or exit the program)
        std::cerr << "Failed to open the wirecell debug log file!" << std::endl;
        return 1; // Return an error code
    }

    if (simplescan) {
       if (world.rank()==0) fmt::print(stderr, "Start simple scan\n");
       fc_master.foreach([world, _mat_frac, _mat_trans, _mat_dirt, _mat_data, dataspec, dim1, dim2, dim3, f_out, rankwork, tol, factr, iter, pearson, xcoord, ycoord, zcoord, &debugFile, debug ](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doScan(b, cp, world.rank(), _mat_frac, _mat_trans, _mat_dirt, _mat_data, dataspec, dim1, dim2, dim3, f_out, rankwork, tol, factr, iter, pearson, xcoord, ycoord, zcoord, debugFile, debug); });
    }/*
    else {
       if (world.rank()==0) fmt::print(stderr, "Start FC -- new version\n");
       fc_master.foreach([world, _mat_frac, _mat_trans, _mat_dirt, _mat_data, dim1, dim2, dim3, nUniverses, rankwork, tol, factr, iter, degree, debug, outer_gpt, inner_gpt, nowrite, msg_every](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doFC(b, cp, world.rank(), _mat_frac, _mat_trans, _mat_dirt, _mat_data, dim1, dim2, dim3 f_out, rankwork, nUniverses, outer_gpt, inner_gpt, tol, factr, iter, degree, debug, doFullGrid, nowrite, msg_every); });
    }

    double endtime   = MPI_Wtime();
    double T13   = MPI_Wtime();
    float _FCtime = endtime-starttime;
    float _modeltime = T11-T10;
    float _tottime = T13-T0;
    float FCtime_max, FCtime_min, FCtime_sum;
    float modeltime_max, tottime_max;
    MPI_Reduce(&_FCtime, &FCtime_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_FCtime, &FCtime_min, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_FCtime, &FCtime_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_FCtime, &FCtime_min, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_modeltime, &modeltime_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_tottime, &tottime_max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);


    if(mfa) std::cout << "This is the MFA implementation!" << std::endl; 
    if (world.rank()==0) fmt::print(stderr, "[{}] that took {} ... {} seconds. Total FC time: {} seconds \n", world.rank(), FCtime_min, FCtime_max, endtime-starttime);
    if (world.rank()==0) fmt::print(stderr, "[{}] building signal model took : {} seconds \n", world.rank(), modeltime_max);
    if (world.rank()==0) fmt::print(stderr, "[{}] Total Time: {} seconds \n", world.rank(), tottime_max);
    if (world.rank()==0) fmt::print(stderr, "Output written to {}\n",out_file);
    */
    //std::cout << "f" << std::endl;
    debugFile.close();
    delete f_out;

    return 0;
}
