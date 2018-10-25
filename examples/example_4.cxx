/*
Example 4: Oscillation Studies
October 19, 2018

Davio Cianci
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "prob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
	std::string xml = "example.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	bool numudis = false;
	bool combined = false;
	int mass_start = -1;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"gen",	no_argument, 0, 'g'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:g", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'g':
				gen = true;
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}
	}

	std::string tag = "EXAMPLE4";

	// This study has two parts. When you first run the code, you must add the -g tag to generate the mass
	// oscillation library. Afterwards, run it without that tag to calculate the chi2.
	// NOTE: This example uses a fullosc method for oscillations, so your xml and input files must have
	// a fullosc channel containing the numu sample fully convolved with nue flux and cross sections.

	
	// PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
	if(gen){
		
		// Load up central value spectra from example 1 and remove the fullosc--this will serve as our background later
        // FOR THIS EXAMPLE TREATING LEESIGNAL AS FULLOSC!!
		SBNspec bkg("EXAMPLE1.SBNspec.root",xml);
		bkg.Scale("leesignal",0.0);
		bkg.WriteOut(tag+"_Bkg");

		// Model has parameters mnu, ue4, um4
		// Since we're precomputing the amplitudes only, we don't really care about the u's
		double mnu = 1.7; //MeV
		NeutrinoModel testModel(mnu,1,1);
	
		// On construction, this will make 3 SBNspecs: 1 sin amplitude, 1  sin2 amplitude and 1 CV oscillated
		SBNgenerate * gen = new SBNgenerate(xml,testModel);

		// Write them to  file
		gen->WritePrecomputedOscSpecs(tag);
	}

	else{
		
		// Load our unoscillated background
		SBNspec bkg(tag+"_Bkg.SBNspec.root",xml);

		// Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
		TFile * fsys = new TFile("EXAMPLE1.SBNcovar.root","read");
		TMatrixD * cov = (TMatrixD*)fsys->Get("frac_covariance_EXAMPLE1");

		// Create our chi2 object with this covariance matrix
		SBNchi my_chi2(bkg,*cov);		

		// Create our desired oscillation model.
		// NOTE: since mass amplitudes are generated above, make sure mnu is same, else we will look for an amplitude
		// that doesn't exist and that won't work.
		NeutrinoModel testModel(1.7, .05, .05);
		SBNosc osc(tag+"_Bkg.SBNspec.root",xml);
		osc.LoadModel(testModel);
		osc.OscillateThis(tag);

		double chi2 = my_chi2.CalcChi(&osc);
		std::cout << "chi2: " << chi2 << std::endl;
	}
	
	return 0;
}
