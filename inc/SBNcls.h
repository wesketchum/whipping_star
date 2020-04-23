#ifndef SBNCLS_H_
#define SBNCLS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNchi.h"
#include "SBNconfig.h"

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

#include "ngrid.h"
#include <gsl/gsl_randist.h>


namespace sbn{


class SBNcls{

	public:

	SBNcls(SBNspec *inh0, SBNspec * inh1, TMatrixD matin) : h0(inh0), h1(inh1), covariance_matrix(matin), chi_h0(*inh0, matin),chi_h1(*inh1,matin){
		which_sample = 0; //default Poisson
        which_mode = 1; //default delta chi
        use_CNP=false;
        maxchival = 210;
		rangen= new TRandom3(0);
        draw_pseudo_from_collapsed = false;
        m_tolerance = 1e-12;
	}
	SBNcls(SBNspec *inh0, SBNspec * inh1) : h0(inh0), h1(inh1), chi_h0(*inh0),chi_h1(*inh1){
		which_sample = 0; //default Poisson
        which_mode = 1; //default delta chi
        use_CNP = false;
        maxchival = 210;
		rangen= new TRandom3(0);
        draw_pseudo_from_collapsed = false;
        m_tolerance = 1e-12;
	}



	SBNspec * h0;
	SBNspec * h1;
	
	SBNchi chi_h0;//previously just chi
	SBNchi chi_h1;

	TMatrixD covariance_matrix;
    bool m_tolerance;

	TRandom3 * rangen;

    bool use_CNP;
    int which_mode;
	int which_sample;
    double maxchival;
	bool draw_pseudo_from_collapsed;
        /****************** Member Functions *************/
    int SetTolerance(double epsilon){
        m_tolerance = epsilon;
        chi_h0.m_tolerance = m_tolerance;            
        chi_h1.m_tolerance = m_tolerance;            
    };
    int SetSampleFromCollapsed(){draw_pseudo_from_collapsed = true;};
    int CalcCLS(int,std::string);
	int SetSampleCovariance();
	int SetSamplePoisson();
    double pval2sig(double p);
    double pval2sig1sided(double p);
    double pval2sig2sided(double p);
    int DrawSampleCovariance(std::string);

    int setMode(int);
    int makePlots(CLSresult &h0_result, CLSresult & h1_result, std::string tag,  int which_mode=0);
    int runConstraintTest();


};


};
#endif
