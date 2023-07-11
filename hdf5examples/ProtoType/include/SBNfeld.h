#ifndef SBNFELD_H_
#define SBNFELD_H_

#include <memory>
#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNchi.h"
#include "SBNconfig.h"
#include "SBNgenerate.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TMath.h"
#include "TGraph.h"

#include "TMath.h"
#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

#include "prob.h"
#include "ngrid.h"
#include <gsl/gsl_randist.h>


namespace sbn{


    class SBNfeld: public SBNconfig{


        NGrid m_grid;

        std::vector<std::vector<double>> m_vec_grid;
        std::vector<SBNspec*> m_cv_spec_grid;
        std::vector<SBNchi*> m_sbnchi_grid;

        TMatrixT<double>  m_full_fractional_covariance_matrix;

        std::unique_ptr<SBNosc>            m_core_spectrum;
        std::unique_ptr<SBNosc>            m_background_spectrum;
        std::unique_ptr<SBNchi>            m_background_chi;
        std::unique_ptr<TVectorT<double> > m_tvec_background_spectrum;

//        bool m_bool_core_spectrum_set;
        bool m_bool_stat_only;
        bool m_bool_print_comparasons;

        int m_num_universes;
        double m_random_seed;
        std::string tag;


        public:
        int nGridPoints() const { return m_grid.f_num_total_points;}
        
        SBNfeld(SBNfeld const&) = delete;
        SBNfeld& operator=(SBNfeld const&) = delete;
        
        SBNfeld(SBNfeld &&) = delete;
        SBNfeld& operator=(SBNfeld &&) = delete;

        virtual ~SBNfeld() {
          //delete m_full_fractional_covariance_matrix;

          for (auto x : m_cv_spec_grid)    delete x;
          for (auto x : m_sbnchi_grid)  delete x;

        } 

        SBNfeld(NGrid const & ingrid, std::string const & intag,  std::string const & inxmlname) : SBNconfig(inxmlname), m_grid(ingrid), tag(intag) {
            m_vec_grid = m_grid.GetGrid();
            m_bool_stat_only = false;
            m_num_universes = 2500;
            m_bool_print_comparasons = false;;
            m_random_seed = -1;
        }
        
        SBNfeld(NGrid const & ingrid, std::string const & intag,  const char * xmldata) : SBNconfig(xmldata, true), m_grid(ingrid), tag(intag) {
            m_vec_grid = m_grid.GetGrid();
            m_bool_stat_only = false;
            m_num_universes = 2500;
            m_bool_print_comparasons = false;;
            m_random_seed = -1;
        }


        //Member Functions
        int FullFeldmanCousins();
        int PointFeldmanCousins(size_t);
        int GlobalScan();
        int RasterScan(); 

        bool statOnly() { return m_bool_stat_only; }
        double seed() {return m_random_seed; }
        TMatrixT<double>   fullFracCovMat() { return m_full_fractional_covariance_matrix; }
        TVectorT<double> *  bgSpectrum() { return m_tvec_background_spectrum.get(); }
        int bgBinsCompressed() { return m_background_spectrum->num_bins_total_compressed; } 


        int GenerateOscillatedSpectra();
        int LoadPreOscillatedSpectra();
        //std::vector<double> LoadPreOscillatedSpectrum(size_t ipoint);
        std::unique_ptr<SBNspec> LoadPreOscillatedSpectrum(size_t ipoint);
        std::unique_ptr<SBNspec> LoadPreOscillatedSpectrum(size_t ipoint, const char * xmldata);
        std::unique_ptr<SBNspec> LoadPreOscillatedSpectrum(size_t ipoint, const char * xmldata,
            std::unordered_map <std::string, std::vector<TH1D> > const & sinsqmap,
            std::unordered_map <std::string, std::vector<TH1D> > const & sinmap, bool compressed=false);

        int SetRandomSeed(double);

        int GenerateBackgroundSpectrum(); 
        int LoadBackgroundSpectrum();
        int LoadBackgroundSpectrum(const char *);
        int LoadBackgroundSpectrum(const char *, std::vector<TH1D> const &);

        int CalcSBNchis();

        int SetCoreSpectrum(std::string const &);
        int SetCoreSpectrum(std::string const &, const char *);
        int SetCoreSpectrum(std::vector<TH1D> histbg, const char * xmldata);

        //int SetFractionalCovarianceMatrix(TMatrixT<double> *);
        int SetFractionalCovarianceMatrix(TMatrixT<double> const & in);
        int SetFractionalCovarianceMatrix(std::string, std::string);
        int SetEmptyFractionalCovarianceMatrix();
        int SetNumUniverses(int);
        int SetStatOnly();

        NeutrinoModel convert3p1(std::vector<double> const & ingrid);


        //This is a stopgap for better SBNchi integration.Hrump, need to fix that wierd float oddity. 
        float CalcChi(std::vector<float>& data, std::vector<double>& prediction, TMatrixT<double> & inverse_covariance_matrix );

    };




}
#endif
