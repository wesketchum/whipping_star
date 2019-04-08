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

#include <gsl/gsl_randist.h>


namespace sbn{


class SBNcls{

	public:

	SBNcls(SBNspec *inh0, SBNspec * inh1, TMatrixD matin) : h0(inh0), h1(inh1), covariance_matrix(matin), chi(*inh0, matin){
		which_sample = 0; //default Poisson
        maxchival = 210;
		rangen= new TRandom3(0);
	}
	SBNcls(SBNspec *inh0, SBNspec * inh1) : h0(inh0), h1(inh1), chi(*inh0){
		which_sample = 0; //default Poisson
        maxchival = 210;
		rangen= new TRandom3(0);
	}



	SBNspec * h0;
	SBNspec * h1;
	
	SBNchi chi;
	TMatrixD covariance_matrix;

	TRandom3 * rangen;

	int which_sample;
    double maxchival;
	/****************** Member Functions *************/
	int CalcCLS(int,std::string);
	int SetSampleCovariance();
	int SetSamplePoisson();
    double pval2sig(double p);
    int DrawSampleCovariance(std::string);
};


};
#endif
