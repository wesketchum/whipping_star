#ifndef SBNCOVARIANCE_H_
#define SBNCOVARIANCE_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include <ctime>

#include "SBNspec.h"
#include "SBNconfig.h"
#include "SBNchi.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TDirectory.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TText.h"
#include "TROOT.h"
#include "TRint.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THnSparse.h"
#include "TStopwatch.h"

#include "params.h"

namespace sbn{

  class SBNcovariance : public SBNconfig {

    bool is_small_negative_eigenvalue;

  public:
		
    SBNspec spec_central_value;	

    std::vector<TH2D> h_spec2d_CV;
    std::vector<TH2D> h_spec2d_sig;
	
    SBNcovariance(std::string xmlname);

    int FormCovarianceMatrix(std::string tag);
    int WriteOut();
    int PrintMatricies(std::string tag);

    int plot_one(TMatrixD matrix, std::string tag, TFile *fin,bool,bool);
    int qualityTesting();
    virtual bool EventSelection(int file);
    virtual int FillHistograms(int file, int uni, double wei);
    
    std::vector<SBNspec> multi_sbnspec;
    std::vector<std::vector<double>> multi_vecspec;

    TMatrixD full_covariance;
    TMatrixD frac_covariance;
    TMatrixD full_correlation;

    std::vector<TMatrixD> vec_full_covariance;
    std::vector<TMatrixD> vec_frac_covariance;
    std::vector<TMatrixD> vec_full_correlation;

    //for plotting, hsa been superseeded by PrintMatricies a bit
    TH2D * hist_frac_cov;
    TH2D * hist_full_cor;
    TH2D * hist_full_cov;

    //Some checks on montecarlos
    double tolerence_positivesemi;
    int universes_used;
    double abnormally_large_weight;

    //Multisim input variables
    std::vector<std::string> variations;
    std::vector<std::vector<double> > vars;

    std::vector<int> num_universes_per_variation;
    std::map<int, std::string> map_universe_to_var;
    std::vector<int> vec_universe_to_var;
    std::map<std::string, int> map_var_to_num_universe;

    std::map<std::string, int> map_var_to_matrix;
    std::vector<int> vec_var_to_matrix;

    int PrintVariations(std::string tag);
    std::string bnbcorrection_str;

    int num_files;

    std::vector<TFile*> files;
    std::vector<TTree*> trees;

    std::vector<int> nentries;
    std::vector<TBranch*>* branch_weight;
    std::vector<std::map<std::string, std::vector<double> >* > f_weights;

    std::vector<std::vector<int> > vars_i;
    std::vector<std::vector<double> > vars_d;

    // In testing 2D fits with 4D covariance..
    SBNspec template_spec;	
    SBNspec spec_central_value2;	
    THnSparseD *frac_covariance_4d;
    TMatrixD full_covariance2D;
    TMatrixD frac_covariance2D;
    TMatrixD full_correlation2D;
    std::vector<std::vector<double> > multi_vecspec2D;
    std::vector<double> vecspec2DCV;
    TStopwatch watch;
    void ProcessEvent(const std::map<std::string, std::vector<double> >& thisfWeight,
		      size_t fileid,
		      int entryid);

  };


};
#endif
