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
#include <regex>
#include <ctime>
#include <chrono>

#include "TMatrixT.h"
#include "TH1D.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNosc.h"
#include "SBNfeld.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
//#include "tools.h"
#include "prob.h"
#include "ngrid.h"

//#undef basic_string_view
using namespace std;
using namespace std::chrono;


struct FitResult {
    size_t n_iter;
    int best_grid_point;
    double last_chi_min, delta_chi;
};

typedef std::array<double, 3> GridPoint;


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


inline Eigen::MatrixXd calcCovarianceMatrixFast(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec) {
    Eigen::MatrixXd ret(M.cols(), M.cols());
    ret.array()    = M.array()*(spec*spec.transpose()).array();
    ret.diagonal() += spec;
    return ret;
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



// Cholesky decomposition and solve for inverted matrix --- fastest
inline Eigen::MatrixXd invertMatrixEigen3(Eigen::MatrixXd const & M){
    return M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
}


inline Eigen::MatrixXd updateInvCov(Eigen::MatrixXd const & covmat, Eigen::VectorXd const & spec_full, sbn::SBNconfig const & conf) {
    auto const & cov = calcCovarianceMatrixFast(covmat, spec_full);
    auto const & out = collapseDetectors(cov, conf);
    return invertMatrixEigen3(out);
}


inline double calcChi(Eigen::VectorXd const & data, Eigen::VectorXd const & prediction, Eigen::MatrixXd const & C_inv ) {
    auto const & diff = data-prediction;
    return diff.transpose() * C_inv * diff;
}
inline double calcChi(Eigen::VectorXd const & diff, Eigen::MatrixXd const & C_inv ) {
    return diff.transpose() * C_inv * diff;
}



class ScaleGenerator {
    public:
        ScaleGenerator(
                sbn::SBNconfig const & conf, std::vector<std::vector<double>> const & vec_grid, Eigen::VectorXd const & core, std::vector<int> const & scale_indicies) 
        {
            m_conf = conf;
            m_gridpoints=vec_grid;
            m_core = core;
            retVec = Eigen::VectorXd(m_conf.num_bins_total);
            m_indicies = scale_indicies;
        }

        Eigen::VectorXd predict(size_t i_grid, bool compressed) {

            //Replace this with a function which scales the predefined signal component. Thats all

            retVec = m_core;
            for(auto &i: m_indicies) retVec[i]*=retVec[i]*m_gridpoints.front()[i_grid];

            if (compressed) return collapseVectorEigen(retVec, m_conf);
            else return retVec;
        }

        size_t gridsize() {return m_gridpoints.front().size();}

    private:
        sbn::SBNconfig m_conf;
        std::vector<std::vector<double>> m_gridpoints;
        Eigen::VectorXd m_core, retVec;
        std::vector<int>m_indicies;
};

inline std::tuple<double, int> universeChi2(Eigen::VectorXd const & data, Eigen::MatrixXd const & C_inv,
        ScaleGenerator signal)
{
    double chimin=std::numeric_limits<double>::infinity();
    Eigen::VectorXd diff(data.rows());
    int bestP(0);
    for (size_t i=0; i<signal.gridsize(); ++i) {
        diff = data - signal.predict(i, true);
        double chi = calcChi(diff, C_inv);
        std::cout<<"   iter "<<i<<" chi: "<<chi<<" sum: "<<diff.sum()<<" "<<bestP<<std::endl;
        if (chi<chimin) {
            chimin = chi;
            bestP=i;
        }
    }
    return {chimin, bestP};
}

inline FitResult performIterativeEigenFit(Eigen::VectorXd const & fake_data, Eigen::VectorXd const & v_coll,
        ScaleGenerator signal,
        Eigen::MatrixXd const & INVCOV,
        Eigen::MatrixXd const & covmat,
        sbn::SBNconfig const & myconf,
        double chi_min_convergance_tolerance = 0.001,
        size_t max_number_iterations = 5
        )
{
    float last_chi_min = FLT_MAX;
    int best_grid_point = -99;
    size_t n_iter = 0;

    Eigen::MatrixXd invcov = INVCOV;//std::vector<double> temp;

    for(n_iter = 0; n_iter < max_number_iterations; n_iter++){
        if(n_iter!=0){
            //Calculate current full covariance matrix, collapse it, then Invert.
            auto const & temp  = signal.predict(best_grid_point, false);
            invcov = updateInvCov(covmat, temp, myconf);
        }
        
        float chi_min = FLT_MAX;
        auto const & resuni  = universeChi2(fake_data, invcov, signal);
        chi_min = std::get<0>(resuni);
        best_grid_point = std::get<1>(resuni);
        if(n_iter!=0){
            if(fabs(chi_min-last_chi_min)< chi_min_convergance_tolerance){
                last_chi_min = chi_min;
                break;
            }
        }
        last_chi_min = chi_min;
    } // End loop over iterations

    //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
    float this_chi = calcChi(fake_data, v_coll, invcov);

    FitResult fr = {n_iter, best_grid_point, last_chi_min, this_chi-last_chi_min}; 
    return fr;
}

Eigen::MatrixXd calcMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
    Eigen::MatrixXd ret(M.cols(), M.cols());
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


Eigen::VectorXd sample(Eigen::VectorXd const & spec, Eigen::MatrixXd const & LMAT, std::mt19937 & rng) {

    std::normal_distribution<double> dist_normal(0,1);
    Eigen::VectorXd RDM(spec.rows());
    for (int i=0;i<spec.rows();++i) RDM[i] = dist_normal(rng);

    return LMAT*RDM + spec;
}

Eigen::VectorXd poisson_fluctuate(Eigen::VectorXd const & spec, std::mt19937 & rng) {
    Eigen::VectorXd RDM(spec.rows());
    for (int i=0;i<spec.rows();++i) {
        std::poisson_distribution<int> dist_pois(spec[i]);
        RDM[i] = double(dist_pois(rng));
    }
    return RDM;
}


void doFC(const char * xmldata, sbn::SBNconfig const & myconf, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
        ScaleGenerator signal, int i_grid,  int nUniverses,  double tol, size_t iter, bool debug, TFile *f, int msg_every=100)
{

    f->cd();

    double tree_delta_chi = 0;
    double tree_chi_min = 0;
    double tree_bf_value=0;
    int tree_bf_pt=0;

    TTree t_outtree(("ttree_"+std::to_string(i_grid)).c_str(),("ttree_"+std::to_string(i_grid)).c_str());
    t_outtree.Branch("delta_chi2",&tree_delta_chi);
    t_outtree.Branch("chi2_min",&tree_chi_min);
    t_outtree.Branch("bf_gridvalue",&tree_bf_value);       
    t_outtree.Branch("bf_gridpoint",&tree_bf_pt);       


    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi, v_nevents;

    results.reserve(nUniverses);
    v_grid.reserve(nUniverses);
    v_univ.reserve(nUniverses);
    v_iter.reserve(nUniverses);
    v_best.reserve(nUniverses);
    v_last.reserve(nUniverses);
    v_dchi.reserve(nUniverses);
    v_nevents.reserve(nUniverses);

    system_clock::time_point t_init = system_clock::now();

    auto const & specfull_e = signal.predict(i_grid, false);
    auto const & speccoll = collapseVectorEigen(specfull_e, myconf);

    std::mt19937 rng(0); // Mersenne twister
    Eigen::MatrixXd const & LMAT = cholD(ECOV, specfull_e);

    for (int uu=0; uu<nUniverses;++uu) {
        std::cout<<"Uni "<<uu<<std::endl;
        auto const & fake_data = poisson_fluctuate(sample(specfull_e, LMAT, rng), rng);//
        auto const & fake_dataC = collapseVectorEigen(fake_data, myconf); 
        std::cout<<fake_data.sum()<<" "<<fake_dataC.sum()<<" "<<fake_data.size()<<" "<<fake_dataC.size()<<std::endl;
        
        results.push_back(performIterativeEigenFit(fake_dataC, speccoll, signal, INVCOVBG, ECOV, myconf, tol, iter));

        v_univ.push_back(uu);
        v_grid.push_back(i_grid);
        v_nevents.push_back(fake_data.sum());

        tree_delta_chi = results.back().delta_chi;
        tree_chi_min = results.back().last_chi_min;
        tree_bf_pt = (int)results.back().best_grid_point;
        t_outtree.Fill();

    }
    system_clock::time_point now = system_clock::now();

    auto t_elapsed = now - t_init;
    auto t_togo = t_elapsed * (1.0 - i_grid)/(i_grid+1);
    //auto t_eta = now + t_togo;
    //std::time_t t = system_clock::to_time_t(t_eta);

    for (auto res : results) {
        v_iter.push_back(res.n_iter);
        v_best.push_back(res.best_grid_point);
        v_last.push_back(res.last_chi_min);
        v_dchi.push_back(res.delta_chi);
    }



    f->cd();
    t_outtree.Write();



}


std::vector<double> flattenHistos(std::vector<TH1D> const & v_hist) {
    std::vector<double> ret; 
    for (auto h : v_hist) {
        for (int i=1; i<(h.GetSize()-1); ++i) ret.push_back(h.GetBinContent(i));
    }
    return ret;
}
void releaseVec(std::vector<double> & vec) {
    vec.clear();
    vec.shrink_to_fit();
}

TMatrixD readFracCovMat(std::string const & rootfile){
    TFile  fsys(rootfile.c_str(),"read");
    TMatrixD cov =  *(TMatrixD*)fsys.Get("frac_covariance");
    fsys.Close();
    std::cout<<"COV: "<<cov.GetNcols()<<" "<<cov.GetNrows()<<std::endl;

    for(int i =0; i<cov.GetNcols(); i++) {
        for(int j =0; j<cov.GetNrows(); j++) {
            if ( std::isnan( cov(i,j) ) )  cov(i,j) = 0.0;
        }
    }
    return cov;
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


int sbn::rootEigenWrapper(std::string const & xmlname, std::string const & tag, std::string const & signal_subchannel, double bkg_scale_value, NGrid & mygrid, int number , bool bool_stat_only){

    double m_chi_min_convergance_tolerance = 0.001;
    sbn::SBNconfig const myconf(xmlname, false);

    std::cout<<1<<std::endl;

    //Core spectrum
    sbn::SBNosc sbncore(tag+"_CV.SBNspec.root",xmlname);
    sbncore.CalcFullVector();
    std::vector<double> core = sbncore.full_vector;
    Eigen::Map<Eigen::VectorXd> ecore(core.data(), core.size(), 1);

    //What incicies should be scaled?
    std::vector<int> scale_indicies = sbncore.GetIndiciesFromSubchannel(signal_subchannel);

    std::cout<<2<<std::endl;

    //BKG spectrum
    sbncore.Scale(signal_subchannel, bkg_scale_value);
    sbncore.CalcFullVector();
    std::vector<double> bkg = sbncore.full_vector;
    Eigen::Map<Eigen::VectorXd> ebkg(bkg.data(), bkg.size(), 1);

    std::cout<<3<<std::endl;
    //Next load up covariance matrix
    TMatrixD covmat;
    covmat.ResizeTo(myconf.num_bins_total,myconf.num_bins_total);
    if (!bool_stat_only) {
        TFile  fsys((tag+".SBNcovar.root").c_str(),"read");
        covmat =  *(TMatrixD*)fsys.Get("frac_covariance");
        std::cout<<"COV: "<<covmat.GetNcols()<<" "<<covmat.GetNrows()<<std::endl;

        for(int i =0; i<covmat.GetNcols(); i++) {
            for(int j =0; j<covmat.GetNrows(); j++) {
                if ( std::isnan( covmat(i,j) ) )  covmat(i,j) = 0.0;
            }
        }
        fsys.Close();
    }
    else {
        covmat.Zero();
    }
    std::cout<<core.size()<<" "<<covmat.GetNrows()<<std::endl;
    Eigen::Map<const Eigen::MatrixXd > ECOVMAT(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());

    //Use the BG only inv cov matrix as start point
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bkg);
    Eigen::Map<const Eigen::MatrixXd> ecov(_cov.GetMatrixArray(), _cov.GetNcols(), _cov.GetNcols());
    auto const & _covcol = collapseDetectors(ecov, myconf);
    auto const & INVCOVBG = invertMatrixEigen3(_covcol);

    auto const vecgrid = mygrid.GetGrid();
    //Finally, the signal generator
    std::cout<<"Creating Signal Generator"<<std::endl;
    ScaleGenerator signal(myconf, vecgrid, ecore, scale_indicies);

    TFile *fout =  new TFile(("SBNfeld_Eigen_output_"+tag+".root").c_str(),"recreate");
    std::cout<<"Starting doFC"<<std::endl;
    doFC(xmlname.c_str(), myconf, ECOVMAT, INVCOVBG, signal, 10, number, m_chi_min_convergance_tolerance, 5, false,fout);

    fout->Close();

    return 0;
}
