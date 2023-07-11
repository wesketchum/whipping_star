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
#include <chrono>
#include <ctime>

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
 
#include "TMatrixT.h"
#include "TH1D.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "tools.h"
#include "prob.h"
#include "ngrid.h"

#include <mfa/mfa.hpp>
#include <mfa/block_base.hpp>

//#undef basic_string_view
using namespace std;
namespace fs = std::filesystem;
using namespace std::chrono;

#include "opts.h"

// set input and ouptut precision here, float or double
#if 0
typedef float                          real_t;
#else
typedef double                         real_t;
#endif

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
   
};

typedef std::array<double, 3> GridPoint;

// TODO what if fixed param is exactly 0 --- pow(10, x) is dangerous then
class GridPoints {
public:
  GridPoints() {}
  GridPoints(std::vector<std::vector<double>> const & m_vec_grid, int setZero=-1) {
    for (auto gp : m_vec_grid) {
      GridPoint P = {pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])};
      if (setZero>=0) P[setZero] =0;//_gridpoints.push_back({pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});
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

class SignalGenerator {
public:
  SignalGenerator(
                  sbn::SBNconfig const & conf, std::vector<std::vector<double>> const & vec_grid, 
		  size_t dim1, size_t dim2, size_t dim3,
                  Eigen::MatrixXd  const & sinsq,
                  Eigen::MatrixXd  const & sin,
                  Eigen::VectorXd const & core, int oscmode
                  ) 
  {
    m_conf = conf;
    m_dim2 = dim1;
    m_dim2 = dim2;
    m_dim2 = dim3;
    m_gridpoints=GridPoints(vec_grid);
    m_sinsq = sinsq;
    m_sin = sin;
    m_core = core;
    m_oscmode = oscmode;
    retVec = Eigen::VectorXd(m_conf.num_bins_total);
  }
  
  Eigen::VectorXd predict(size_t i_grid, bool compressed) {
    auto const & gp = m_gridpoints.Get(i_grid);
    sbn::NeutrinoModel this_model(gp[0]*gp[0], gp[1], gp[2], false);
    int m_idx = massindex(i_grid);
    Oscillate(m_sinsq.row(m_idx), m_sin.row(m_idx), this_model);
    
    if (compressed) return collapseVectorEigen(m_core+retVec, m_conf);
      else return m_core+retVec;
  }
  
  int massindex(size_t igrd) {return int(floor( (igrd) / m_dim2 ));}
  size_t gridsize() {return m_gridpoints.NPoints();}
  GridPoint getgrid(size_t igrd) {return m_gridpoints.Get(igrd);};
  
  void Oscillate(Eigen::VectorXd const & sf_sinsq, Eigen::VectorXd const & sf_sin, 
                 sbn::NeutrinoModel & working_model) 
  {
    retVec.setZero(); // !!!
    int which_dm = 41; // FIXME This must not be hardcoded!!! 41 is the value to be used in the 1 sterile case
    
    double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);
    
    // TODO Make this logic nice
    if (m_oscmode==0) {
      prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
      prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
      prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
      prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
    } 
    else if (m_oscmode==1) prob_mumu = working_model.oscAmp(2,2,which_dm,2);
    else if (m_oscmode==2) {
      // This allows for both nu_e dis/app and nu_mu dis
      prob_mumu = working_model.oscAmp(2,2,which_dm,2);
      prob_ee = working_model.oscAmp(1,1,which_dm,2);
      prob_mue = working_model.oscAmp(2,1,which_dm,1);
      prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
      prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
      prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
    }
    else {
      std::cerr << "oscillation mode has to be either 0 or 1: " << m_oscmode << "\n";
      exit(1);
    }
    
    double osc_amp(0), osc_amp_sq(0);
    int osc_pattern(0);
    // Iterate over channels
    size_t offset(0);
    for (int i=0; i<m_conf.num_channels; i++) {
      size_t nbins_chan = m_conf.num_bins[i];
      auto const & thisPattern = m_conf.subchannel_osc_patterns[i];
      for (int j=0; j<m_conf.num_subchannels[i]; j++){
        osc_pattern = thisPattern[j];
        switch (osc_pattern){
        case 11:
          osc_amp_sq = prob_ee;
          osc_amp = 0;
          break;
        case -11:
          osc_amp_sq = prob_ee;
          osc_amp = 0;
          break;
        case 22:
          osc_amp_sq = prob_mumu;
          osc_amp = 0;
          break;
        case -22:
          osc_amp_sq = prob_mumu;
          osc_amp = 0;
          break;
        case 21:
          osc_amp    = prob_mue;
          osc_amp_sq = prob_mue_sq;
          break;
        case -21:
          osc_amp    = prob_muebar;
          osc_amp_sq = prob_muebar_sq;
          break;
        case 0: 
          osc_amp = 0;
          osc_amp_sq = 0;
        default:
          break;
        }
        
        // Iterate over detectors
        for (int d=0; d<m_conf.num_detectors;++d) {
          size_t first  = d*m_conf.num_bins_detector_block + offset;
          retVec.segment(first, nbins_chan).noalias() += osc_amp   *  sf_sin.segment(first, nbins_chan);
          retVec.segment(first, nbins_chan).noalias() += osc_amp_sq*  sf_sinsq.segment(first, nbins_chan);
        }
        offset +=nbins_chan;
      }
    }
  }
  
private:
  sbn::SBNconfig m_conf;
  GridPoints m_gridpoints;
  size_t m_dim1, m_dim2, m_dim3;
  Eigen::MatrixXd m_sinsq, m_sin;
  Eigen::VectorXd m_core, retVec;
  int m_oscmode;
};


TMatrixD readFracCovMat(std::string const & rootfile){
  TFile  fsys(rootfile.c_str(),"read");
  TMatrixD cov =  *(TMatrixD*)fsys.Get("frac_covariance");
  fsys.Close();
  
  for(int i =0; i<cov.GetNcols(); i++) {
    for(int j =0; j<cov.GetNrows(); j++) {
      if ( std::isnan( cov(i,j) ) )  cov(i,j) = 0.0;
    }
  }
  return cov;
}

std::vector<TH1D> readHistos(std::string const & rootfile, std::vector<string> const & fullnames) {
  std::vector<TH1D > hist;
  TFile f(rootfile.c_str(),"read");
  for (auto fn: fullnames){ hist.push_back(*((TH1D*)f.Get(fn.c_str())));}
  f.Close();
  return hist;
}

std::vector<double> flattenHistos(std::vector<TH1D> const & v_hist) {
  std::vector<double> ret; 
  for (auto h : v_hist) {
    for (int i=1; i<(h.GetSize()-1); ++i) ret.push_back(h.GetBinContent(i));
  }
  return ret;
}

std::vector<std::tuple<std::string, float> > getFilesAndDm(std::string const & inDir, std::string const & tag, std::string const & subthing, double const & m_min, double const & m_max, bool debug=false) {
  const std::regex re("[-+]?([0-9]*\\.[0-9]+|[0-9]+)");
  std::smatch match;
  std::string result;
  
  std::vector<std::tuple<std::string, float> > ret;
  
  for (const auto & entry : fs::directory_iterator(inDir)) {
    //std::cout << "I'm here 1: " << inDir << "+" <<tag << "+" << subthing << std::endl;
    if (std::string(fs::path(entry).stem()).rfind(tag + subthing, 0) == 0) {
      //std::cout << "I'm here 2" << std::endl;
      std::string test     = std::string(fs::path(entry));
      std::string teststem = std::string(fs::path(entry).stem());
      
      std::size_t loc = teststem.find(subthing);
      std::string _test = teststem.substr(loc);
      //std::cout << "test, teststem, _test: " << test << ", " << teststem << ", " << _test << std::endl;
      if (std::regex_search(_test, match, re) && match.size() > 1) {
        float lgmsq = std::stof(match.str(0));
        if (lgmsq > m_max || lgmsq < m_min) {
          if (debug) std::cerr << "\t NOT using file " << test << " with " << match.str(0) << " " << lgmsq << "\n";
          continue;
        }
        ret.push_back({test, lgmsq});
        if (debug) std::cerr << "\t Using file " << test << " with " << match.str(0) << " " << lgmsq << "\n";
      }
    }
  }
  return ret;
}

std::tuple< std::vector<std::vector<double>>, std::vector<float>> mkHistoVecStd(std::string const & inDir, std::string const & tag, std::string const & objname, std::vector<string> const & fullnames, double const & m_min, double const & m_max, bool debug=false ) {
  
  // The input files come unordered
  auto const & inputs = getFilesAndDm(inDir, tag, objname, m_min, m_max, debug);
  std::vector<std::vector<double> > temp(inputs.size());
  
  if (inputs.size() ==0) {
    std::cerr << "Error, no valid input files, exiting. Maybe check --xmin --xmax and location -i\n";
    exit(1);
  }
  // Sort by mass, ascending
  std::vector<float> masses;
  for (auto in : inputs) masses.push_back(std::get<1>(in));
  std::sort(masses.begin(), masses.end());
  
  std::cerr << "Summary of mass dimension --- we have " << inputs.size() << " inputs corresponding to these mass-squared splittings: \n";
  for (auto m : masses) std::cerr << m << " ";
  std::cerr << "\n";
  for (auto in : inputs) {
    std::vector<float>::iterator it = std::find(masses.begin(), masses.end(), std::get<1>(in));
    int mass_index = std::distance(masses.begin(), it);
    temp[mass_index] = flattenHistos(readHistos(std::get<0>(in), fullnames));
  }
   
  return {temp, masses};
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

// Cholesky decomposition and solve for inverted matrix --- fastest
inline Eigen::MatrixXd invertMatrixEigen3(Eigen::MatrixXd const & M){
  return M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
}

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

void createDataSets(HighFive::File* file, size_t nPoints, size_t full, size_t coll) {
  file->createDataSet<int>   ("i_grid",       HighFive::DataSpace( { nPoints,       1} ));
  file->createDataSet<double>("delta_msq",    HighFive::DataSpace( { nPoints,       1} ));
  file->createDataSet<double>("Ue4",          HighFive::DataSpace( { nPoints,       1} ));
  file->createDataSet<double>("Umu4",         HighFive::DataSpace( { nPoints,       1} ));
  file->createDataSet<double>("specfull",     HighFive::DataSpace( { nPoints,    full} ));
  file->createDataSet<double>("speccoll",     HighFive::DataSpace( { nPoints,    coll} ));
}

void createH5InputSpectrums(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
                            sbn::SBNconfig const & myconf, Eigen::MatrixXd const & ECOV, 
                            Eigen::MatrixXd const & INVCOVBG, Eigen::VectorXd const & ecore, 
                            SignalGenerator signal, std::vector<size_t> const & rankwork, bool debug)
{
  
  int nUniverses=1;
  double starttime, stoptime, endtime;
  size_t nbinsfull = ecore.size();
  size_t nbinscoll = collapseVectorEigen(ecore, myconf).size();
  
  //not sure if the work balancing is  needed for making the spectrum, 
  //might be necessary for high-dimensionality and fine resolution 
  size_t pStart = rankwork[0];
  size_t uStart = rankwork[1];
  size_t pLast  = rankwork[2];
  size_t uLast  = rankwork[3];
  
  size_t i_begin = pStart * uStart;
  size_t i_end   = pLast  * uLast;
  
  //fmt::print(stderr, "[{}] a,b,c,d: {} {} {} {} start at {} end at {}  lends {}\n", rank, pStart, uStart, pLast, uLast, i_begin, i_end, i_end-i_begin);
  size_t lenDS = i_end - i_begin;
  
  Eigen::VectorXi mat_grid(lenDS);
  Eigen::VectorXd mat_delta_msq(lenDS); 
  Eigen::VectorXd mat_Ue4(lenDS); 
  Eigen::VectorXd mat_Umu4(lenDS);
  Eigen::MatrixXd mat_specfull(lenDS,nbinsfull); 
  Eigen::MatrixXd mat_speccoll(lenDS,nbinscoll);
  
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
  
  starttime = MPI_Wtime();
  
  for (auto r : RW) {
    int i_grid = r[0];
    
    auto const & gp       = signal.getgrid(i_grid);
    auto const & specfull = signal.predict(i_grid,false);
    auto const & speccoll = collapseVectorEigen(specfull, myconf);
    
    //Eigen Matrices to store the spectrums and relative oscillation parameters
    //example for 3N+1 oscillation
    mat_grid(i_grid) = i_grid;
    mat_delta_msq(i_grid) = gp[0]*gp[0];
    mat_Ue4(i_grid) = gp[1];
    mat_Umu4(i_grid) = gp[2];
    mat_specfull.row(i_grid) = specfull;
    mat_speccoll.row(i_grid) = speccoll;
  }
  
  stoptime = MPI_Wtime();
  
  if (rank==0){ 
    fmt::print(stderr, "Total time to create the input spectrums: {}", stoptime-starttime);

    //H5Easy::DumpOptions options(H5Easy::Compression(9), H5Easy::DumpMode::Overwrite);

    H5Easy::File file1("inputfiles_spectrumgridpoint.h5", H5Easy::File::Overwrite);
    H5Easy::dump(file1, "spec_full", mat_specfull);
    H5Easy::dump(file1, "spec_coll", mat_speccoll);
    H5Easy::dump(file1, "i_grid", mat_grid.transpose());
    H5Easy::dump(file1, "delta_msq", mat_delta_msq.transpose());
    H5Easy::dump(file1, "Ue4", mat_Ue4.transpose());
    H5Easy::dump(file1, "Umu4", mat_Umu4.transpose());

    file1.flush();
  }
  endtime = MPI_Wtime();
  if (rank==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", rank, endtime-starttime);
  std::cout << "exiting function" << std::endl;
}

inline bool file_exists (const std::string& name) {
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
  
  //std::cerr << MPI_COMM_WORLD.size() <<"\n";
  //createTimingDataSets(f_time);
  
  size_t nPoints=-1;
  int mode=0;
  int msg_every=100;
  size_t nUniverses=1;
  int NTEST(0);
  std::string out_file="inputfilespectrums.f5";
  std::string f_BG="NuMuDis_BKG_ONLY.SBNspec.root";
  std::string f_CV="NuMuDis_CV.SBNspec.root";
  std::string f_COV="NuMuDis.SBNcovar.root";
  std::string tag="";
  std::string d_in="";
  std::string xml="";
  double xmin(-1.0);
  double xmax(1.1);
  double xwidth(0.1);
  double ymin(-2.3);
  double ymax(0.1);
  double ywidth(0.05);
  double zmin(-2.3);
  double zmax(0.1);
  double zwidth(0.05);
  
  // get command line arguments
  using namespace opts;
  Options ops(argc, argv);
  ops >> Option('x', "xml",        xml,        "XML config.");
  ops >> Option('t', "tag",        tag,        "Tag.");
  ops >> Option('i', "indir",      d_in,       "Input file directory.");
  ops >> Option("core",            f_CV,       "Central values filename.");
  ops >> Option('b', "background", f_BG,       "Backgrounds filename.");
  ops >> Option('c', "covmat",     f_COV,      "Covariance matrix filename.");
  ops >> Option("xmin",            xmin,       "xmin");
  ops >> Option("xmax",            xmax,       "xmax");
  ops >> Option("xwidth",          xwidth,     "xwidth");
  ops >> Option("ymin",            ymin,       "ymin");
  ops >> Option("ymax",            ymax,       "ymax");
  ops >> Option("ywidth",          ywidth,     "ywidth");
  ops >> Option("zmin",            zmin,       "zmin");
  ops >> Option("zmax",            zmax,       "zmax");
  ops >> Option("zwidth",          zwidth,     "zwidth");
  ops >> Option("mode",            mode,       "Mode 0 is default, mode 1 is muon disappearance mode 2 is muon disapp/electron app");
  ops >> Option('o', "output",     out_file,   "Output filename."); 
  bool debug       = ops >> Present('d', "debug", "Operate on single gridpoint only");
  
  if( world.rank()==0 ) {
    fmt::print(stderr, "\n*** This is diy running SBN Feldman Cousins ***\n");
  }
  
  // Whole bunch of tests
  if ( world.rank() == 0 ) {
    std::cout << "f_BG, f_COV, f_CV, xml: " << f_BG << ", " << f_COV << ", " << f_CV << ", " << xml << std::endl;
    std::vector<std::string> infiles = {f_BG, f_COV, f_CV, xml};
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
    if (d_in=="") {
      std::cerr << "Input dir (-i, --indir) cannot be undefined, exiting\n";
      exit(1);
    }
  }
  
  double T1   = MPI_Wtime();
  double time0 = MPI_Wtime();
  
  // Read the xml file on rank 0
  std::string line, text;
  if ( world.rank() == 0 ) {
    std::ifstream in(xml);
    while(std::getline(in, line))  text += line + "\n";
  }
  // YUCK, is this really the most elegant way to broadcast a simple string???
  int textsize = text.size();
  MPI_Bcast(&textsize, 1, MPI_INT, 0, world);
  if ( world.rank() != 0 ) text.resize(textsize);
  MPI_Bcast(const_cast<char*>(text.data()), textsize, MPI_CHAR, 0, world);
    
  double T2   = MPI_Wtime();
  // Central configuration object
  const char* xmldata = text.c_str();
  sbn::SBNconfig myconf(xmldata, false);
  
  double T3   = MPI_Wtime();
  // Pre-oscillated spectra
  std::vector<double> sinsqvec, sinvec;
  std::vector<float> msqsplittings;
  std::vector<Eigen::VectorXd > sinsqvec_eig, sinvec_eig;
  int nFilesIn(0);
  if (world.rank()==0) {
    auto temp = mkHistoVecStd(d_in, tag, "_SINSQ_", myconf.fullnames, xmin, xmax, debug);
    std::cout << "****" << std::endl;
    sinsqvec = asVector(std::get<0>(temp));
    msqsplittings = std::get<1>(temp);
    
    auto temp2 = mkHistoVecStd(d_in, tag, "_SIN_", myconf.fullnames, xmin, xmax, debug);
    sinvec   = asVector(std::get<0>(temp2));
    if (sinsqvec.size() != sinvec.size()) {
      std::cerr << "Error, number of input files for _SINSQ_ (" << sinsqvec.size() << ") differs from _SIN_ (" << sinvec.size() << ") exiting.\n";
      exit(1);
    }
    nFilesIn = msqsplittings.size();
    
  }
  double T4   = MPI_Wtime();
  diy::mpi::broadcast(world, sinsqvec, 0);
  diy::mpi::broadcast(world, sinvec,   0);
  diy::mpi::broadcast(world, msqsplittings,   0);
  diy::mpi::broadcast(world, nFilesIn, 0);
  
  double T5   = MPI_Wtime();
  for (auto v : splitVector(sinsqvec, nFilesIn)) sinsqvec_eig.push_back(Eigen::Map<Eigen::VectorXd> (v.data(), v.size(), 1) );
  for (auto v : splitVector(sinvec, nFilesIn)) sinvec_eig.push_back(Eigen::Map<Eigen::VectorXd> (v.data(), v.size(), 1) );
  Eigen::Map<Eigen::VectorXf> masses(msqsplittings.data(), msqsplittings.size(), 1);
  
  int nrows, ncols, ncols_coll;
  nrows = nFilesIn;
  ncols = sinsqvec_eig[0].size();
  ncols_coll = collapseVectorEigen(sinsqvec_eig[0],myconf).size();
  Eigen::MatrixXd sinsqeig(nrows,ncols);
  Eigen::MatrixXd sineig(nrows,ncols);
  Eigen::MatrixXd sinsqeigcoll(nrows,ncols_coll);
  Eigen::MatrixXd sineigcoll(nrows,ncols_coll);
  
  for (int i = 0; i < sinsqvec_eig.size(); ++i){
    sinsqeig.row(i) = Eigen::VectorXd::Map(&sinsqvec_eig[i][0], sinsqvec_eig[0].size());
    sinsqeigcoll.row(i) = collapseVectorEigen(sinsqeig.row(i), myconf);
    sineig.row(i) = Eigen::VectorXd::Map(&sinvec_eig[i][0], sinvec_eig[0].size());
    sineigcoll.row(i) = collapseVectorEigen(sineig.row(i), myconf);
  }
  
  // Core spectrum
  std::cout << "f_CV, f_BG, f_COV: " << f_CV << ", " << f_BG << ", " << f_COV << std::endl;
  std::vector<double> core;
  if (world.rank()==0) {
    auto const & cvhist = readHistos(f_CV, myconf.fullnames);
    core = flattenHistos(cvhist);
  }
  diy::mpi::broadcast(world, core, 0);
  
  Eigen::Map<Eigen::VectorXd> ecore(core.data(), core.size(), 1);
  // Background
  std::vector<double> bgvec;
  if (world.rank() == 0) {
    std::vector<TH1D> bghist  = readHistos(f_BG, myconf.fullnames);
    bgvec = flattenHistos(bghist);
    for (int i=0; i<myconf.num_channels; i++) {
      size_t nbins_chan = myconf.num_bins[i];
      for (int j=0; j<myconf.num_subchannels[i]; j++){
        if(myconf.subchannel_names[i][j]=="fullosc"){
          for (int k=0; k<nbins_chan; k++){
            int a = i*nbins_chan+j*nbins_chan+k;
            bgvec[a]=0.0;
          }
        }
      }
    }
    bghist.clear();
    bghist.shrink_to_fit();
  }
  diy::mpi::broadcast(world, bgvec, 0);
  Eigen::Map<Eigen::VectorXd> ebg(bgvec.data(), bgvec.size(), 1);
  
  // Read the covariance matrix on rank 0 --- broadcast and subsequently buid from array
    TMatrixD covmat;
    
    size_t nBins(0);
    std::vector<double> v_covmat;
    
    if ( world.rank() == 0 ) {
      TMatrixD temp = readFracCovMat(f_COV);
      nBins = temp.GetNcols();
      const double *pData = temp.GetMatrixArray();
      v_covmat.assign(pData, pData + temp.GetNoElements());
    }
    // broadcast
    diy::mpi::broadcast(world, v_covmat, 0);
    diy::mpi::broadcast(world, nBins,    0);
    // Set data of TMatrix
    covmat.ResizeTo(nBins, nBins);
    covmat.SetMatrixArray(v_covmat.data());
    releaseVec(v_covmat);
    
    Eigen::Map<const Eigen::MatrixXd > ECOVMAT(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());
    
    // Use the BG only inv cov matrix as start point
    // Use the core inv cov matrix as start point or sensitivity
    //TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    std::cout << "core: " << ecore << std::endl;
    std::cout << "bgvec: "; 
    for (int b=0; b < bgvec.size(); b++ ) std::cout << bgvec[b] << std::endl;
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    Eigen::Map<const Eigen::MatrixXd> ecov(_cov.GetMatrixArray(), _cov.GetNcols(), _cov.GetNcols());
    auto const & _covcol = collapseDetectors(ecov, myconf);
    auto const & INVCOVBG = invertMatrixEigen3(_covcol);
    
    //Now let's write the grid for the input spectrums
    
    // Fetch the masses (or the first dimension) and count the number of files for the grid resolution
    
    size_t dim1(0);
    size_t dim2(0);
    size_t dim3(0);
    size_t nDims(0);
    
    // Setup grid
    NGrid mygrid;
    
    double mmin = masses.head(1)(0)*0.5;
    double mmax = masses.tail(1)(0)*0.5;
    double mwidth = (masses(1) - masses(0))*0.5;
    
    if (world.rank()==0) std::cerr << "Mass setup for input grid: " << mmin << " < " << mmax << " width: " << mwidth << "\n";
    
    int setZero=-1;
    
    mygrid.AddDimension("m4", mmin, mmax, mwidth );
    dim1=mygrid.f_dimensions[0].f_N;
    if (mode==0) {
      mygrid.AddDimension("ue4", ymin, ymax, ywidth);// arbirtrarily dense! mixing angle nu_e
      mygrid.AddFixedDimension("um4", 0.0);
      setZero=2;
      dim2 = mygrid.f_dimensions[1].f_N;
      dim3 = 1;
      nDims = 2;
    }
    else if (mode==1) {
      mygrid.AddFixedDimension("ue4", 0.0);
      mygrid.AddDimension("um4", ymin, ymax, ywidth);
      setZero=1;
      dim2 = mygrid.f_dimensions[2].f_N;
      dim3 = 1;
      nDims = 2;
    }
    else if (mode==2) {
      mygrid.AddDimension("ue4", ymin, ymax, ywidth);
      mygrid.AddDimension("um4", zmin, zmax, zwidth);
      dim2 = mygrid.f_dimensions[1].f_N;
      dim3 = mygrid.f_dimensions[2].f_N;
      nDims = 3;
    }
    else {
      std::cerr << "Error, the mode must be either 0 or 1 a the moment: " << mode << "\n";
      exit(1);
    }
    
    nPoints = mygrid.f_num_total_points;
    GridPoints GP(mygrid.GetGrid(), setZero);
    
    //create the oscillation spectrun OR any other physics model spectrum
    SignalGenerator signal(myconf, mygrid.GetGrid(), dim1, dim2, dim3, sinsqeig, sineig, ebg, mode);
    
    if(world.rank()==0){
      H5Easy::File file("corespectrum.h5", H5Easy::File::Overwrite);
      H5Easy::dump(file, "bgvec", ebg);
      H5Easy::dump(file, "core", ecore);
      H5Easy::dump(file, "covmat", ECOVMAT);
      H5Easy::dump(file, "invcovbg", INVCOVBG);
    }
    if(world.rank()==0){
      std::cout << "make _SINSQ_.h5 " << std::endl;
      H5Easy::File file("_SINSQ_.h5", H5Easy::File::Overwrite);
      H5Easy::dump(file, "_SINSQ_", sinsqeig);
      H5Easy::dump(file, "_SINSQCOLL_", sinsqeigcoll);
      H5Easy::dump(file, "masses", msqsplittings);
    }
    if(world.rank()==0){
      H5Easy::File file("_SIN_.h5", H5Easy::File::Overwrite);
      H5Easy::dump(file, "_SIN_", sineig);
      H5Easy::dump(file, "_SINCOLL_", sineigcoll);
      H5Easy::dump(file, "masses", msqsplittings);
    }
    
    // Setup the blocks for parallelization. 
    // Useful for high dimensional problem with large number of grid points
    // a block refers to a gridpoint
    size_t blocks = world.size();//nPoints;
    if (world.rank()==0) fmt::print(stderr, "creation of the input spectrum will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds<real_t> spectrum_domain(1);
    spectrum_domain.min[0] = 0.;
    spectrum_domain.max[0] = blocks-1;
    
    diy::FileStorage               storage("./DIY.XXXXXX");
    diy::RoundRobinAssigner        spectrum_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds<real_t>> spectrum_decomposer(1, spectrum_domain, blocks);
    diy::RegularBroadcastPartners  spectrum_comm(    spectrum_decomposer, 1, true);
    diy::RegularMergePartners      spectrum_partners(spectrum_decomposer, 1, true);
    diy::Master                    spectrum_master(world, 1, -1, &Block<real_t>::create, &Block<real_t>::destroy, &storage, &Block<real_t>::save, &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    spectrum_decomposer.decompose(world.rank(), assigner,
                                  [&](int gid, const Bounds<real_t>& core, const Bounds<real_t>& bounds, const Bounds<real_t>& domain, const RCLink<real_t>& link)
                                  { Block<real_t>::add(gid, core, bounds, domain, link, spectrum_master, 1, nBins+1, 0.0); });
    
    //do LoadBalance
    size_t _S(nPoints);
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
    
    //write out the input spectrums
    if (world.rank()==0) fmt::print(stderr, "Write input spectrums \n");
    spectrum_master.foreach([world, ECOVMAT, INVCOVBG, ecore, myconf, rankwork, debug, signal ](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    { createH5InputSpectrums(b, cp, world.rank(), myconf, ECOVMAT, INVCOVBG, ecore, signal, rankwork, debug); });
    //createH5InputSpectrums( world.rank(), myconf, ECOVMAT, INVCOVBG, ecore, signal, rankwork, debug);
    return 0;
}
