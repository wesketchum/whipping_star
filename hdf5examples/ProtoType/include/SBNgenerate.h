#ifndef SBNGENERATE_H_
#define SBNGENERATE_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "SBNspec.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "prob.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TROOT.h"
#include "TRint.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THnSparse.h"

#include <map>
#include <ctime>
#include "params.h"

namespace sbn{


class SBNgenerate : public SBNconfig{
	

	public:
		
	SBNspec spec_central_value;	
	NeutrinoModel nu_model;

	SBNspec spec_osc_sin;
	SBNspec spec_osc_sinsq;
	
	SBNgenerate(std::string xmlname, NeutrinoModel inModel );
	SBNgenerate(std::string xmlname);

	int WriteOut();

	virtual bool EventSelection(int file);
	virtual int FillHistograms(int file, int uni, double wei);
	
	int WritePrecomputedOscSpecs(std::string tag);

	//Some checks on multisims

	//Multisim input variables
	std::vector<std::vector<double> > vars;

	int num_files;
	std::vector<TFile *> files;
	std::vector<TTree *> trees;

	std::vector<int> nentries;
	std::vector< TBranch *> * branch_weight;
    std::vector<std::map<std::string, std::vector<double> >* > f_weights;


	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;


};


};
#endif
