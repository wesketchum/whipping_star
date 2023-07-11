#pragma GCC optimize("O3","unroll-loops","inline", "peel-loops")

#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
 
#include "TMatrixT.h"
#include "TH1D.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
//#include "nomad.hpp"


using namespace std;

#include "opts.h"

typedef diy::DiscreteBounds Bounds;

struct FitResult {
   size_t n_iter;
   int best_grid_point;
   double last_chi_min, delta_chi;
};






struct Block
{
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    std::vector<double> last_chi_min, delta_chi;
    std::vector<int> best_grid_point, n_iter, i_grid, i_univ;
    std::vector< std::vector<float> > fake_data;
    std::vector< std::vector<double> > spec_full, spec_coll, inv_cov;
};

// Split input vector into pieces of similar size
std::vector<std::vector<int>> splitVector(const std::vector<int>& vec, size_t n) {
    std::vector<std::vector<int>> outVec;
    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;
    size_t begin = 0;
    size_t end = 0;

    for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
        end += (remain > 0) ? (length + !!(remain--)) : length;
        outVec.push_back(std::vector<int>(vec.begin() + begin, vec.begin() + end));
        begin = end;
    }
    return outVec;
}

// Concatenate vector of vectors as vector
std::vector<double> asVector(std::vector<std::vector<double> > v_in) {
    std::vector<double> allSpectra;
    allSpectra.reserve(v_in.size()*v_in[0].size());
    for (auto temp : v_in) allSpectra.insert(allSpectra.end(), temp.begin(), temp.end());
    return allSpectra;
}

std::vector<double> asVector(TMatrixT<double> const & M) {
    const double *pData = M.GetMatrixArray();
    std::vector<double> vData;
    vData.assign(pData, pData + M.GetNoElements());
    return vData;
}

template<typename T>
std::vector<T> myslice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;
    std::vector<T> vec(first, last);
    return vec;
}

std::vector< std::vector<double> > sliceVector(std::vector<double> const & input, int nPoints) {
   std::vector< std::vector<double> > test;
   test.reserve(nPoints);
   int nBins = input.size()/nPoints;

   std::vector<double> work;
   work.reserve(nBins);

   for (int i=0;i<nPoints;++i) {
      work=myslice(input, i*nBins, (i+1)*nBins-1);
      test.push_back(work);
   }
   return test;
}

// Properly release memory
void releaseVec(std::vector<double> & vec) {
    vec.clear();
    vec.shrink_to_fit();
}

// Reduce loaded spectra and inverted matrices, also write to hdf5
void gatherSpectra(Block* b, const diy::ReduceProxy& rp, const diy::RegularMergePartners& partners,
         std::vector<std::vector<double> > & full_out,
         std::vector<std::vector<double> > & coll_out,
         std::vector<std::vector<double> > & invcov_out,
         HighFive::File* file,
         bool otf) {
    // step 1: dequeue and merge
    for (int i = 0; i < rp.in_link().size(); ++i) {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid()) continue;
        std::vector< std::vector<double> > in_vals_spec_full, in_vals_spec_coll, in_vals_inv_cov;
        rp.dequeue(nbr_gid, in_vals_spec_full);
        rp.dequeue(nbr_gid, in_vals_spec_coll);
        if (!otf) rp.dequeue(nbr_gid, in_vals_inv_cov);
        for (size_t j=0; j<in_vals_spec_full.size(); ++j) (b->spec_full).push_back(in_vals_spec_full[j]);
        for (size_t j=0; j<in_vals_spec_coll.size(); ++j) (b->spec_coll).push_back(in_vals_spec_coll[j]);
        if (!otf) {for (size_t j=0; j<in_vals_inv_cov.size();   ++j) (b->inv_cov  ).push_back(in_vals_inv_cov[j]);}
    }

    // step 2: enqueue
    for (int i = 0; i < rp.out_link().size(); ++i) {
        // only send to root of group, but not self
        if (rp.out_link().target(i).gid != rp.gid()) rp.enqueue(rp.out_link().target(i), b->spec_full);
        if (rp.out_link().target(i).gid != rp.gid()) rp.enqueue(rp.out_link().target(i), b->spec_coll);
        if (!otf) {if (rp.out_link().target(i).gid != rp.gid()) rp.enqueue(rp.out_link().target(i), b->inv_cov);}
    }

    // step 3: write out
    if (rp.gid()==0 && rp.out_link().size()==0) {
       HighFive::DataSet spec = file->getDataSet("specfull");
       spec.select(         {0, 0}, {size_t(b->spec_full.size()), size_t(b->spec_full[0].size())}).write(b->spec_full);
       HighFive::DataSet coll = file->getDataSet("speccoll");
       coll.select(         {0, 0}, {size_t(b->spec_coll.size()), size_t(b->spec_coll[0].size())}).write(b->spec_coll);

       full_out.reserve(b->spec_full.size());
       for (auto ff : b->spec_full) full_out.push_back(ff);
       coll_out.reserve(b->spec_coll.size());
       for (auto ff : b->spec_coll) coll_out.push_back(ff);
       if (!otf) {
          HighFive::DataSet invc = file->getDataSet("invcov");
          invc.select(         {0, 0}, {size_t(b->inv_cov.size()),   size_t(b->inv_cov[0].size())}).write(b->inv_cov);
          invcov_out.reserve(b->inv_cov.size());
          for (auto ff : b->inv_cov) invcov_out.push_back(ff);
       }

    }
}

void createDataSets(HighFive::File* file, size_t nPoints, size_t nBins, size_t nBinsC, size_t nUniverses) {
   //HighFive::DataSetCreateProps props;
   //props.add(HighFive::Deflate(4));
    file->createDataSet<double>("specfull",     HighFive::DataSpace( { nPoints           , nBins  } ));
    file->createDataSet<double>("speccoll",     HighFive::DataSpace( { nPoints           , nBinsC } ));
    file->createDataSet<double>("invcov",       HighFive::DataSpace( { nPoints+1         , nBinsC*nBinsC } ));
    file->createDataSet<double>("fakedata",     HighFive::DataSpace( { nPoints*nUniverses, nBinsC } ));
    file->createDataSet<double>("last_chi_min", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<double>("delta_chi",    HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("best_grid_point", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("n_iter",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    // Some bookkeeping why not
    file->createDataSet<int>("i_grid",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<int>("i_univ",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<double>("gridx",        HighFive::DataSpace( {nPoints,                   1} ));
    file->createDataSet<double>("gridy",        HighFive::DataSpace( {nPoints,                   1} ));
}

void writeGrid(HighFive::File* file, std::vector<std::vector<double> > const & coords) {
    std::vector<double> xcoord;
    std::vector<double> ycoord;

    for (size_t i=0; i< coords.size(); i++) {
       xcoord.push_back(coords[i][0]);
       ycoord.push_back(coords[i][1]);
    }
    HighFive::DataSet d_gridx          = file->getDataSet("gridx");
    HighFive::DataSet d_gridy          = file->getDataSet("gridy");
    d_gridx.select(   {0, 0}, {xcoord.size(), 1}).write(xcoord);
    d_gridy.select(   {0, 0}, {ycoord.size(), 1}).write(ycoord);
}

void writeInvCovBG(HighFive::File* file, TMatrixT<double> const & M, size_t pos) {
    std::vector<double> vInvCovBG = asVector(M);
    HighFive::DataSet invc = file->getDataSet("invcov");
    invc.select(     {size_t(pos), 0}, {1, size_t(vInvCovBG.size())}).write(vInvCovBG);
}

//
TMatrixD readFracCovMat(std::string const & rootfile){
    TFile  fsys(rootfile.c_str(),"read");
    TMatrixD cov =  *(TMatrixD*)fsys.Get("frac_covariance");
    fsys.Close();
    return cov;
}


std::vector<TH1D> readHistos(std::string const & rootfile, std::vector<string> const & fullnames) {
    std::vector<TH1D > hist;
    TFile f(rootfile.c_str(),"read");
    for (auto fn: fullnames) hist.push_back(*((TH1D*)f.Get(fn.c_str())));
    f.Close();
    return hist;
}

std::unordered_map <std::string, std::vector<TH1D> > mkHistoMap(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::unordered_map <std::string, std::vector<TH1D> > hmap;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      hmap[mass_tag] = readHistos(fname, fullnames);
   }
   return hmap;
}

// Assume order is correct
std::vector<double> flattenHistos(std::vector<TH1D> const & v_hist) {
   std::vector<double> temp;
   for (auto hist : v_hist) {
      for (int i=1; i< (hist.GetNbinsX()+1); ++i) {
         temp.push_back(hist.GetBinContent(i));
      }
   }
   return temp;
}
//std::vector<TH1D> vectToHist(std::vector<double> const & flat) {
   //std::vector<TH1D> ret;
   //ret.resize(flat.size());
   //for (auto vec : flat) {
      //TH1D temp();
      //ret.push_back(temp):
   //}
//}

float calcChi(std::vector<float> const & data, std::vector<double> const & prediction, TMatrixT<double> const & C_inv ) {
    int n = data.size();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > CINV(C_inv.GetMatrixArray(), n, n);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > PRED(prediction.data(), n, 1);
    // Massage data right now to convert from float to double.
    std::vector<double> temp;
    temp.assign(data.begin(), data.end());
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > DATA(temp.data(), n, 1);
    Eigen::VectorXd DIFF = DATA-PRED;
    return DIFF.transpose() * CINV * DIFF; 
}

inline size_t grid_to_index(size_t ix, size_t iy) {
   // Total grid size static for now ( number of y-values in this 2D case)
   size_t NY(48);
   return ix*NY + iy;
}


//class My_Evaluator : public NOMAD::Evaluator {
//public:
  //My_Evaluator  ( const NOMAD::Parameters & p, std::vector<float> const & data, std::vector<std::vector<double> > const & predictions, TMatrixT<double> const & C_inv  ) :
    //NOMAD::Evaluator ( p ) {_data = data; _predictions=predictions;_C_inv.ResizeTo(C_inv.GetNcols(), C_inv.GetNcols());_C_inv=C_inv;}

  //~My_Evaluator ( void ) {}

  //bool eval_x ( NOMAD::Eval_Point   & x          ,
		//const NOMAD::Double & h_max      ,
		//bool                & count_eval   ) const
  //{

     //int ix = x[0].round();
     //int iy = x[1].round();
     ////std::cerr << x[0] << " " << x[1] << " becomes " << ix << " " << iy << " ###\n";
     //float val = calcChi(_data, _predictions[grid_to_index(ix, iy)], _C_inv);
     ////std::cerr << x[0] << " " << x[1] << " becomes " << ix << " " << iy <<  " with val " << val << " ###\n";

    //x.set_bb_output  ( 0 , val  ); // objective value

    //count_eval = true; // count a black-box evaluation

    //return true;       // the evaluation succeeded
  //}

  //std::vector<float>  _data;
  //std::vector<std::vector<double> >  _predictions;
  //TMatrixT<double>  _C_inv;
	
	
//};

// Calculate all chi2 for this universe
std::vector<double> universeChi2(std::vector<float> const & data, std::vector<std::vector<double> > const & predictions, TMatrixT<double> const & C_inv ) {
   std::vector<double> result;
   //std::vector<double> result2;
   result.reserve(predictions.size());
   //result2.reserve(predictions.size());

   //double starttime, endtime;
   //starttime = MPI_Wtime();
   for (auto p : predictions) {
      result.push_back(calcChi(data, p, C_inv));
   }

   return result;
}
   //float chi_min = *std::min_element(result.begin(), result.end());
   //int  best_grid_point = std::min_element(result.begin(), result.end()) - result.begin();
   //std::cerr << "BEST GLOBAL: " << chi_min << " at " << best_grid_point << "\n";
   //endtime   = MPI_Wtime(); 
   //fmt::print(stderr, "standard took {} seconds\n", endtime-starttime);
   //NOMAD::Display out ( std::cout );
   //out.precision ( NOMAD::DISPLAY_PRECISION_STD );
   ////NOMAD::begin (0, wtf);

   //// parameters creation:
   //NOMAD::Parameters p ( out );

   //p.set_DIMENSION (2);             // number of variables

   //vector<NOMAD::bb_output_type> bbot (3); // definition of
   //bbot[0] = NOMAD::OBJ;                   // output types
   //bbot[1] = NOMAD::PB;
   //bbot[2] = NOMAD::EB;
   //p.set_BB_OUTPUT_TYPE ( bbot );

////   p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
   //p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
   //p.set_BB_INPUT_TYPE ( 0 , NOMAD::INTEGER );
   //p.set_BB_INPUT_TYPE ( 1 , NOMAD::INTEGER );
   //p.set_LOWER_BOUND ( NOMAD::Point ( 2 , 0 ) );
   //NOMAD::Point ub ( 2 );
   //ub[0]=20;
   //ub[1]=47;

   //p.set_UPPER_BOUND ( ub );

   //p.set_X0 ( NOMAD::Point(2 , 5) );  // starting point

   //p.set_MAX_BB_EVAL (100);     // the algorithm terminates after
                                //// 100 black-box evaluations
   //p.set_DISPLAY_DEGREE(2);
   ////p.set_SOLUTION_FILE("sol.txt");

   //// parameters validation:
   //p.check();

   //// custom evaluator creation:
   //My_Evaluator ev   ( p,  data, predictions, C_inv);

   //starttime = MPI_Wtime();
   //// algorithm creation and execution:
   //NOMAD::Mads mads ( p , &ev );
   //mads.run();
   //endtime   = MPI_Wtime(); 
   //const NOMAD::Eval_Point * sol = mads.get_best_feasible();
   //int sol_x = sol->value(0);
   //int sol_y = sol->value(1);
   ////delete sol;
   //fmt::print(stderr, "MADS took {} seconds: {}  {}  {}  {}\n", endtime-starttime, sol_x, sol_y, grid_to_index(sol_x, sol_y), sol->get_f());


   ////exit(1);
   ////for (size_t ix=0;ix<21;++ix) {
      ////for (size_t iy=0;iy<48;++iy) result2.push_back(calcChi(data, predictions[grid_to_index(ix, iy)], C_inv));
   ////}

   ////for (size_t i=0;i<result.size(); ++i) std::cerr << result[i] << " vs " << result2[i] << " : " << result[i] - result2[i] << "\n";


   //return result;
//}

TMatrixT<double> calcCovarianceMatrix(TMatrixT<double> const & M, std::vector<double> const & spec){
    TMatrixT<double> Mout( M.GetNcols(), M.GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M.GetNcols(); i++) {
        for(int j =0; j<M.GetNrows(); j++) {
            if ( std::isnan( M(i,j) ) )  Mout(i,j) = 0.0;
            else                         Mout(i,j) = M(i,j)*spec[i]*spec[j];

            if (i==j) Mout(i,i) += spec[i];
        }
    }
    return Mout;
}

// Only used for covariance matrix a very beginning, just being lazy here
TMatrixT<double> invertMatrix(TMatrixT<double> &M){
    TMatrixT<double> McI(M.GetNrows(),M.GetNrows());
    McI.Zero();
    TDecompSVD svd(M);
    if (!svd.Decompose()) {
        std::cerr<<" (InvertMatrix) Decomposition FAILED, matrix not symettric?, has nans?" << std::endl;
        std::cerr<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< M.GetNrows(); i++){
            for(int j=0; j< M.GetNrows(); j++) std::cerr<<i<<" "<<j<<" "<<M(i,j)<<std::endl;
        }
    }

    else McI = svd.Invert();

    if (!McI.IsValid()) std::cout << "ERROR: The inverted matrix isnt valid! Something went wrong.." << std::endl;
    return McI;
}

// SVD --- slowest
Eigen::MatrixXd invertMatrixEigen(TMatrixT<double> &M, double epsilon = std::numeric_limits<double>::epsilon()){
    int n = M.GetNrows();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > COV(M.GetMatrixArray(), n, n);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(COV,Eigen::ComputeThinU|Eigen::ComputeThinV );
    double tolerance = epsilon * std::max(COV.cols(), COV.rows()) *svd.singularValues().array().abs()(0);
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

// LU decomposition
Eigen::MatrixXd invertMatrixEigen2(TMatrixT<double> &M){
    int n = M.GetNrows();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > COV(M.GetMatrixArray(), n, n);
    return COV.inverse();
}

// Cholesky decomposition and solve for inverted matrix --- fastest
Eigen::MatrixXd invertMatrixEigen3(TMatrixT<double> &M){
    int n = M.GetNrows();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > COV(M.GetMatrixArray(), n, n);
    return COV.llt().solve(Eigen::MatrixXd::Identity(n,n));
}


// This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void collapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){

    std::vector< std::vector< TMatrixT<double> > > Summed(conf.num_channels, std::vector<TMatrixT<double>>(conf.num_channels) );	//Initialise a matrix of matricies, to ZERO.
    for(int ic = 0; ic < conf.num_channels; ic++){
        for(int jc =0; jc < conf.num_channels; jc++){
            Summed[ic][jc].ResizeTo(conf.num_bins[jc], conf.num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
            Summed[ic][jc] = 0.0;
        }
    }

    int mrow = 0.0;
    int mcol = 0.0;

    for(int ic = 0; ic < conf.num_channels; ic++){ 	 // Loop over all rows
        for(int jc =0; jc < conf.num_channels; jc++){    // Loop over all columns
            for(int m=0; m < conf.num_subchannels[ic]; m++){
                for(int n=0; n< conf.num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing together
                   int a, b, c, d;
                   a=mrow + n*conf.num_bins[jc];
                   b=mrow + n*conf.num_bins[jc] + conf.num_bins[jc]-1;
                   c=mcol + m*conf.num_bins[ic];
                   d=mcol + m*conf.num_bins[ic] + conf.num_bins[ic]-1;
                   Summed[ic][jc] +=  M.GetSub(a,b,c,d);
                }
            }
            mrow += conf.num_subchannels[jc]*conf.num_bins[jc]; // As we work our way left in columns, add on that many bins
        } // end of column loop
        mrow = 0; // as we end this row, reSet row count, but jump down 1 column
        mcol += conf.num_subchannels[ic]*conf.num_bins[ic];
    } // end of row loop

    ///********************************* And put them back toGether! ************************//
    Mc.Zero();
    mrow = 0;
    mcol = 0;

    //Repeat again for Contracted matrix
    for(int ic = 0; ic < conf.num_channels; ic++){
        for(int jc =0; jc < conf.num_channels; jc++){
            Mc.SetSub(mrow, mcol, Summed[ic][jc]);
            mrow += conf.num_bins[jc];
        }
        mrow = 0;
        mcol +=conf.num_bins[ic];
    }
}

void collapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){
    Mc.Zero();
    int nrow = conf.num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
    int crow = conf.num_bins_detector_block_compressed; //N_e_bins+N_m_bins;
    for (int m=0; m<conf.num_detectors; m++) {
        for (int n=0; n<conf.num_detectors; n++) {
            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);
            imat = M.GetSub(n*nrow, n*nrow+nrow-1, m*nrow, m*nrow+nrow-1);
            collapseSubchannels(imat, imatc, conf);
            Mc.SetSub(n*crow, m*crow, imatc);
        }
    }
}

void collapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){
    Mc.Zero();
    int nrow = conf.num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
    int crow = conf.num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

    for (int m=0; m<conf.num_modes ; m++) {
        for (int n=0; n<conf.num_modes; n++) {
            TMatrixT<double> imat(nrow, nrow);
            TMatrixT<double> imatc(crow, crow);
            imat = M.GetSub(n*nrow, n*nrow+nrow-1, m*nrow, m*nrow+nrow-1);
            collapseDetectors(imat, imatc, conf);
            Mc.SetSub(n*crow, m*crow, imatc);
        }
    }
}

//unique_ptr<sbn::SBNspec> loadPreOscillatedSpectrum(int i_uni, string const & tag, string const & xml, NGrid const & mygrid, const char * xmldata, std::vector<TH1D> const & bghist,  std::vector<TH1D> const & cvhist, TMatrixD const & covmat) {
unique_ptr<sbn::SBNspec> loadPreOscillatedSpectrum(int i_grid, string const & tag, NGrid const & mygrid, const char * xmldata, std::vector<TH1D> const & cvhist) {
    sbn::SBNfeld myfeld(mygrid, tag, xmldata);
    myfeld.SetCoreSpectrum(cvhist, xmldata);
    myfeld.SetRandomSeed(-1);
    return myfeld.LoadPreOscillatedSpectrum(i_grid, xmldata);
};

unique_ptr<sbn::SBNspec> loadPreOscillatedSpectrum(
      int i_grid, string const & tag, NGrid const & mygrid, const char * xmldata,
      std::vector<TH1D> const & cvhist,
      std::unordered_map <std::string, std::vector<TH1D> > const & sinsqmap,
      std::unordered_map <std::string, std::vector<TH1D> > const & sinmap) {

    sbn::SBNfeld myfeld(mygrid, tag, xmldata);
    myfeld.SetCoreSpectrum(cvhist, xmldata);
    myfeld.SetRandomSeed(-1);
    return myfeld.LoadPreOscillatedSpectrum(i_grid, xmldata, sinsqmap, sinmap);
};

//void loadSpectrum(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid const & mygrid, HighFive::File* f_out, const char * xmldata, std::vector<TH1D> const & bghist,  std::vector<TH1D> const & cvhist,TMatrixD const & covmat) {
void loadSpectrum(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank,
      string tag, NGrid const & mygrid, HighFive::File* f_out, const char * xmldata,
      std::vector<TH1D> const & cvhist,
      std::unordered_map <std::string, std::vector<TH1D> > const & sinsqmap,
      std::unordered_map <std::string, std::vector<TH1D> > const & sinmap) {

    auto myspec = loadPreOscillatedSpectrum(cp.gid(), tag, mygrid, xmldata, cvhist, sinsqmap, sinmap);
    std::vector<double> full = myspec->full_vector;
    std::vector<double> coll = myspec->collapsed_vector;
    
    double sum = std::accumulate(full.begin(), full.end(), 0.0);
    if (sum >1e10 || sum <10) {
       std::cerr << "[" << cp.gid() <<"] sum of specfull: " << sum << "\n";
       myspec->PrintFullVector();
       std::cerr << "and the collapsed vector \n";
       myspec->PrintCollapsedVector();
       abort();
    }
    b->spec_full.push_back(full);
    b->spec_coll.push_back(coll);
}
void loadSpectrum(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, NGrid const & mygrid, HighFive::File* f_out, const char * xmldata,  std::vector<TH1D> const & cvhist) {
    auto myspec = loadPreOscillatedSpectrum(cp.gid(), tag, mygrid, xmldata, cvhist);
    std::vector<double> full = myspec->full_vector;
    std::vector<double> coll = myspec->collapsed_vector;
    
    double sum = std::accumulate(full.begin(), full.end(), 0.0);
    if (sum >1e10 || sum <10) {
       std::cerr << "[" << cp.gid() <<"] sum of specfull: " << sum << "\n";
       myspec->PrintFullVector();
       std::cerr << "and the collapsed vector \n";
       myspec->PrintCollapsedVector();
       abort();
    }
    b->spec_full.push_back(full);
    b->spec_coll.push_back(coll);
}

void updateInvCov(TMatrixT<double> & invcov, TMatrixD const & covmat, std::vector<double> const & spec_full, sbn::SBNconfig const & conf, int nBinsC) {
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, spec_full);
    TMatrixT<double> _covcol;
    _covcol.ResizeTo(nBinsC, nBinsC);
    collapseModes(_cov, _covcol, conf);
    Eigen::MatrixXd temp = invertMatrixEigen3(_covcol);
    invcov.SetMatrixArray(temp.data());
}

// DIY 
void invertMatrix(Block* b, diy::Master::ProxyWithLink const& cp, TMatrixD const & covmat, int nBins, int nBinsC, sbn::SBNconfig const & conf) {

    TMatrixT<double> _cov;
    TMatrixT<double> _covcol;
    TMatrixT<double> _covcolinv;
    _cov.ResizeTo(nBins, nBins);
    _covcol.ResizeTo(nBinsC, nBinsC);
    _covcolinv.ResizeTo(nBinsC, nBinsC);

    // At this stage, each block should only ever have one spectrum
    assert(b->spec_full.size()==1);
    _cov = calcCovarianceMatrix(covmat, b->spec_full[0]);
    collapseModes(_cov, _covcol, conf);
    Eigen::MatrixXd temp = invertMatrixEigen3(_covcol);
    _covcolinv.SetMatrixArray(temp.data());

    const double *pData = _covcolinv.GetMatrixArray();
    std::vector<double> v_covmat;
    v_covmat.assign(pData, pData + _covcolinv.GetNoElements());
    b->inv_cov.push_back(v_covmat);
}


FitResult coreFC(std::vector<float> const & fake_data, std::vector<double> const & v_coll,
      std::vector< std::vector<double> > const & spectra,
      std::vector< std::vector<double> > const & allColl,
      std::vector<TMatrixT<double> > const & INVCOV, TMatrixT<double> const & covmat, 
      sbn::SBNconfig const & myconf, int nBinsC, bool onthefly, 
      double chi_min_convergance_tolerance = 0.001,
      size_t max_number_iterations = 5)
{
   float last_chi_min = FLT_MAX;
   int best_grid_point = -99;
   size_t n_iter = 0;
   TMatrixT<double> inv_cov = INVCOV.back();// inverse_background_collapsed_covariance_matrix;
   
   for(n_iter = 0; n_iter < max_number_iterations; n_iter++){
       if(n_iter!=0){
           //Calculate current full covariance matrix, collapse it, then Invert.
           if (onthefly) updateInvCov(inv_cov, covmat, spectra[best_grid_point], myconf, nBinsC);
           else inv_cov = INVCOV[best_grid_point];
       }
       //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
       float chi_min = FLT_MAX;
       std::vector<double> allchi = universeChi2(fake_data, allColl,  inv_cov);
       chi_min = *std::min_element(allchi.begin(), allchi.end());
       best_grid_point = std::min_element(allchi.begin(), allchi.end()) - allchi.begin();
 
       if(n_iter!=0){
           //Step 3.0 Check to see if min_chi for this particular fake_data  has converged sufficiently
           if(fabs(chi_min-last_chi_min)< chi_min_convergance_tolerance){
               last_chi_min = chi_min;
               break;
           }
       }
       last_chi_min = chi_min;
   } // End loop over iterations

   //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
   // float this_chi = this->CalcChi(fake_data, true_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix);
   float this_chi = calcChi(fake_data, v_coll, inv_cov);

   FitResult fr = {n_iter, best_grid_point, last_chi_min, this_chi-last_chi_min}; 
   return fr;
}

void doFC(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, 
      NGrid mygrid, bool dry, bool onthefly, TMatrix const & covmat, const char * xmldata, 
      std::vector<TMatrixT<double> > const & INVCOV, 
      std::vector< std::vector<double> > const & spectra, 
      std::vector< std::vector<double> > const & allColl, 
      sbn::SBNconfig const & myconf, int nUniverses, 
      HighFive::File* file, std::vector<int> const & rankwork,
      double tol, size_t iter, bool debug)
{

    double starttime, endtime;
    std::vector<float> fake_data;
    std::vector<std::vector<float> > FD;
    FD.reserve(rankwork.size()*nUniverses);
    std::vector<FitResult> results;
    results.reserve(rankwork.size()*nUniverses);
   
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    v_grid.reserve(rankwork.size()*nUniverses);
    v_univ.reserve(rankwork.size()*nUniverses);
    v_iter.reserve(rankwork.size()*nUniverses);
    v_best.reserve(rankwork.size()*nUniverses);
    std::vector<double> v_last, v_dchi;
    v_last.reserve(rankwork.size()*nUniverses);
    v_dchi.reserve(rankwork.size()*nUniverses);
    
    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;

       std::vector<double> specfull = spectra[i_grid];
       sbn::SBNspec myspec(specfull, xmldata, i_grid, false);
       
       int nBinsColl = allColl[i_grid].size();

       // NOTE: the true spectrum collapse is allColl[i_grid] --- no need to call myspec.Collapse as SampleCov uses the full spectrum
    
       sbn::SBNchi mychi(myspec, covmat, xmldata, false);//, myfeld.seed()); FIXME SEEED?

       starttime = MPI_Wtime();
       for (int uu=0; uu<nUniverses;++uu) {
          // OK this is indeed stupid, the 
          //double t1 = MPI_Wtime();
          fake_data = mychi.SampleCovariance(&myspec); // YUCK this thing rebuilds the covariance matrix*spectrum all the time
          //double t2 = MPI_Wtime();
          results.push_back(coreFC(fake_data, allColl[i_grid], spectra, allColl, INVCOV, covmat, myconf, nBinsColl, onthefly, tol, iter));
          //double t3 = MPI_Wtime();
          FD.push_back(fake_data);
          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
          //double t4 = MPI_Wtime();
          //if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} sampling took {} seconds\n",cp.gid(), i_grid, rankwork.size(), t2-t1);
          //if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} coreFC took {} seconds\n",cp.gid(), i_grid, rankwork.size(), t3-t2);
          //if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} rest took {} seconds\n",cp.gid(), i_grid, rankwork.size(), t4-t3);
       }
       endtime   = MPI_Wtime(); 
       if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} ({} universes) took {} seconds\n",cp.gid(), i_grid, rankwork.size(), nUniverses, endtime-starttime);
       myspec.reset();
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
    // This is for the fake data dump
    HighFive::DataSet d_fakedata        = file->getDataSet("fakedata");

    size_t d_bgn = rankwork[0]*nUniverses;
    for (auto res : results) {
       v_iter.push_back(res.n_iter);
       v_best.push_back(res.best_grid_point);
       v_last.push_back(res.last_chi_min);
       v_dchi.push_back(res.delta_chi);
    }

    d_fakedata.select(         {d_bgn, 0}, {size_t(FD.size()), size_t(FD[0].size())}).write(FD);
    d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
    d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
    d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
    d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
    d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
    d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
    endtime   = MPI_Wtime();
    if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
}

inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// --- main program ---//
int main(int argc, char* argv[]) {
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBinsC=90;
    int nPoints=-1;
    int mockFactor=1;
    int nUniverses=1;
    std::string out_file="test.hdf5";
    std::string f_BG="BOTHv2_BKG_ONLY.SBNspec.root";
    std::string f_CV="BOTHv2_CV.SBNspec.root";
    std::string f_COV="BOTHv2.SBNcovar.root";
    std::string tag="";
    std::string xml="";
    double xmin(-1.0);
    double xmax(1.1);
    double xwidth(0.1);
    double ymin(-2.3);
    double ymax(0.1);
    double ywidth(0.05);
    double tol(0.001);
    size_t iter(5);
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    ops >> Option('o', "output",     out_file,   "Output filename.");
    ops >> Option('u', "nuniverse",  nUniverses, "Number of universes");
    ops >> Option('x', "xml",        xml,        "XML config.");
    ops >> Option("tol",             tol,        "Minimiser tolerance");
    ops >> Option("iter",            iter,       "Max number of iterations.");
    ops >> Option('t', "tag",        tag,        "Tag.");
    ops >> Option("nbinsC",          nBinsC,     "Number of collapsed bins in 2d dataset");
    ops >> Option("core",            f_CV,       "Central values filename.");
    ops >> Option('b', "background", f_BG,       "Backgrounds filename.");
    ops >> Option('c', "covmat",     f_COV,      "Covariance matrix filename.");
    ops >> Option("xmin",            xmin,       "xmin");
    ops >> Option("xmax",            xmax,       "xmax");
    ops >> Option("xwidth",          xwidth,     "xwidth");
    ops >> Option("ymin",            ymin,       "ymin");
    ops >> Option("ymax",            ymax,       "ymax");
    ops >> Option("ywidth",          ywidth,     "ywidth");
    ops >> Option("mock",            mockFactor, "Mockfactor");
    bool debug       = ops >> Present('d', "debug", "Operate on single gridpoint only");
    bool dryrun      = ops >> Present("dry", "dry run");
    bool onthefly    = ops >> Present("otf", "on the fly matrix inversion (saves memory)");
    bool statonly    = ops >> Present("stat", "Statistical errors only");
    //bool storeINVCOV = ops >> Present("storeinvcov", "Write inverted covariance matrices to hdf5");
    if (ops >> Present('h', "help", "Show help")) {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }
    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
    }

    NGrid mygrid;
    mygrid.AddDimension("m4",  xmin, xmax, xwidth);//0.1 --- can't change easily --- sterile mass
    mygrid.AddDimension("ue4", ymin, ymax, ywidth);//0.1 --- arbirtrarily dense! mixing angle nu_e --- can change ywidth. Try 20^6
    mygrid.AddFixedDimension("um4", 0.0); //0.05

    nPoints = mygrid.f_num_total_points;


    // Whole bunch of tests
    if ( world.rank() == 0 ) {
       if (world.size() > nPoints) {
          std::cerr << "Impossible to run on more ranks than grid points, exiting.\n";
          exit(1);
       }
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
    }

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

    const char* xmldata = text.c_str();
    sbn::SBNconfig myconf(xmldata, false);

    // TODO find a way to serialize and reconstruct TH1D that is not too painful
    double time0 = MPI_Wtime();
    std::unordered_map <std::string, std::vector<TH1D> > sinsqmap = mkHistoMap(tag+"_SINSQ_dm_", myconf.fullnames);
    std::unordered_map <std::string, std::vector<TH1D> > sinmap   = mkHistoMap(tag+"_SIN_dm_",   myconf.fullnames);
    double time1 = MPI_Wtime();
    if (world.rank()==0) fmt::print(stderr, "[{}] loading spectra the new way took {} seconds\n",world.rank(), time1 - time0);

    // Background
    std::vector<double> bgvec;
    if (world.rank() == 0) {
       std::vector<TH1D> bghist  = readHistos(f_BG, myconf.fullnames);
       bgvec = flattenHistos(bghist);
       bghist.clear();
       bghist.shrink_to_fit();
    }
    diy::mpi::broadcast(world, bgvec, 0);

    // Central values
    // FIXME is this done on every rank?
    std::vector<TH1D> cvhist = readHistos(f_CV, myconf.fullnames);

    // Read the covariance matrix on rank 0 --- broadcast and subsequently buid from array
    TMatrixD covmat;

    size_t nBins(0);
    if (!statonly) {
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
    }
    else {
       nBins=bgvec.size();
       covmat.ResizeTo(nBins, nBins);
       covmat.Zero();
    }

    if( world.rank()==0 ) {
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    nBins:  {}\n", nBins);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    f_BG:  {}\n", f_BG);
      fmt::print(stderr, "\n    f_CV:  {}\n", f_CV);
      fmt::print(stderr, "\n    f_COV:  {}\n", f_COV);
      fmt::print(stderr, "\n    iter :  {}\n", iter);
      fmt::print(stderr, "\n    tol:    {}\n", tol);
      if (statonly) fmt::print(stderr, "\n    S T A T  O N L Y \n");
      if (debug) fmt::print(stderr, "\n    D E B U G \n");
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
      fmt::print(stderr, "***********************************\n");
    }
    
    // Create hdf5 file structure here 
    HighFive::File* f_out  = new HighFive::File(out_file,
                        HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));

    // Create datasets needed TODO nbinsC --- can we get that from somewhere?
    createDataSets(f_out, nPoints, nBins, nBinsC, nUniverses);
   
    // First rank also writes the grid so we know what the poins actually are
    if (world.rank() == 0)  writeGrid(f_out, mygrid.GetGrid());

    // First diy setup to load spectra
    size_t blocks = nPoints;
    Bounds domain(1);
    domain.min[0] = 0;
    domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> decomposer(1, domain, blocks);
    diy::RegularBroadcastPartners  comm(    decomposer, 2, true);
    diy::RegularMergePartners      partners(decomposer, 2, true);
    diy::Master master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), domain, assigner, master);

    double starttime, endtime;
    starttime = MPI_Wtime();
    //master.foreach([world, tag, xml, mygrid, f_out, xmldata, bghist, cvhist, covmat](Block* b, const diy::Master::ProxyWithLink& cp)
                           //{loadSpectrum(b, cp, world.size(), world.rank(), tag, xml, mygrid, f_out, xmldata, bghist, cvhist, covmat); });
    //master.foreach([world, tag, mygrid, f_out, xmldata, cvhist](Block* b, const diy::Master::ProxyWithLink& cp)
                           //{loadSpectrum(b, cp, world.size(), world.rank(), tag, mygrid, f_out, xmldata, cvhist); });
    master.foreach([world, tag, mygrid, f_out, xmldata, cvhist, sinsqmap, sinmap](Block* b, const diy::Master::ProxyWithLink& cp)
                           {loadSpectrum(b, cp, world.size(), world.rank(), tag, mygrid, f_out, xmldata, cvhist, sinsqmap, sinmap); });
    endtime = MPI_Wtime();
    if (world.rank()==0) fmt::print(stderr, "[{}] loading spectra took {} seconds\n",world.rank(), endtime-starttime);

    if (!onthefly) {
       starttime = MPI_Wtime();
       master.foreach([covmat, nBins, nBinsC, myconf](Block* b, const diy::Master::ProxyWithLink& cp)
                              {invertMatrix(b, cp, covmat, nBins, nBinsC, myconf); });
       endtime = MPI_Wtime();
       if (world.rank()==0) fmt::print(stderr, "[{}] matrix inversions took {} seconds\n",world.rank(), endtime-starttime);
    }

    std::vector< std::vector<double> > specFull, specColl, invCov;
    diy::reduce(master, assigner, partners,
        [&specFull, &specColl, &invCov, world, f_out, onthefly](Block* b, const diy::ReduceProxy& rp, const diy::RegularMergePartners& partners){
        gatherSpectra(b, rp, partners, specFull, specColl, invCov, f_out, onthefly); });

    // TODO is there merit in deleting the blocks after the reduction???

    std::vector<double> flatFull, flatColl, flatCov;
    if (world.rank()==0) flatFull = asVector(specFull);
    if (world.rank()==0) flatColl = asVector(specColl);

    diy::mpi::broadcast(world, flatFull, 0);
    diy::mpi::broadcast(world, flatColl, 0);

    specFull = sliceVector(flatFull, nPoints);
    specColl = sliceVector(flatColl, nPoints);
    releaseVec(flatFull);
    releaseVec(flatColl);
    
    if (!onthefly) {
       if (world.rank()==0) flatCov  = asVector(invCov);
       diy::mpi::broadcast(world, flatCov, 0);
       invCov   = sliceVector(flatCov , nPoints);
       releaseVec(flatCov);
    }

    // Reconstruct inverted covariance TMatrices from vector double
    starttime = MPI_Wtime();
    std::vector<TMatrixT<double> > INVCOV;

    if (!onthefly) {
       INVCOV.reserve(nPoints*mockFactor);
       for (int ii = 0;ii<mockFactor; ++ii) {
          for (auto ic : invCov) {
             TMatrixT<double> temp;
             temp.ResizeTo(nBinsC, nBinsC);
             temp.SetMatrixArray(ic.data());
             INVCOV.push_back(temp);
          }
       }
    }

    //std::vector<double> bla;
    //for (auto x : bgvec) bla.push_back(1.0);
    //TMatrixT<double> _temp = calcCovarianceMatrix(covmat,bla);
    //TMatrixT<double> _tempcol;
    //_tempcol.ResizeTo(nBinsC, nBinsC);
    //starttime = MPI_Wtime();
    //for (int i=0;i<1000;++i) {
       //collapseModes(_temp, _tempcol, myconf);
    //}
    //endtime = MPI_Wtime();
    //if (world.rank()==0) fmt::print(stderr, "[{}] 1000 collapse operations took {} seconds\n",world.rank(), endtime-starttime);
    //if (world.rank()==0) {
       //TFile* f = new TFile("bla.root", "RECREATE");
       //_temp.Write("annoyingcrap");
       //f->Close();
       //delete f;
    //}



    
    // Add the BG
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    TMatrixT<double> _covcol;
    _covcol.ResizeTo(nBinsC, nBinsC);
    collapseModes(_cov, _covcol, myconf);
    INVCOV.push_back(invertMatrix(_covcol));
    //std::vector<double> vInvCovBG = asVector(INVCOV.back());
    writeInvCovBG(f_out, INVCOV.back(), size_t(nPoints));
    //HighFive::DataSet invcc = f_out->getDataSet("invcov");
    //invc.select(     {size_t(nPoints), 0}, {1, size_t(vInvCovBG.size())}).write(vInvCovBG);
    endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] matrix broadcast etc took {} seconds\n",world.rank(), endtime-starttime);

    //// Now more blocks as we have universes
    blocks = world.size();//nPoints;//*nUniverses;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds fc_domain(1);
    fc_domain.min[0] = 0;
    fc_domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> fc_decomposer(1, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, 2, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, 2, true);
    diy::Master                    fc_master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), fc_domain, fc_assigner, fc_master);//, share_face, wrap, ghosts);


    std::vector<int> work(nPoints);
    std::iota(std::begin(work), std::end(work), 0); //0 is the starting number
    std::vector<int> rankwork = splitVector(work, world.size())[world.rank()];

    starttime = MPI_Wtime();
    fc_master.foreach([world,  mygrid, dryrun, onthefly, covmat, xmldata, INVCOV, specFull, specColl, myconf, nUniverses, f_out, rankwork, tol, iter, debug](Block* b, const diy::Master::ProxyWithLink& cp)
                           { doFC(b, cp, world.size(), world.rank(), mygrid, dryrun, onthefly, covmat, xmldata, INVCOV, specFull, specColl, myconf, nUniverses, f_out, rankwork, tol, iter, debug); });
    endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] FC took {} seconds\n",world.rank(), endtime-starttime);

    starttime = MPI_Wtime();

    delete f_out;
    return 0;
}
