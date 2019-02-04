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
#include "SBNcls.h"
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
  bool sample_from_covariance = true;
  int num_MC_events = 100000;
  std::string tag = "TEST";

  std::string signal_file = "EMPTY";
  std::string background_file = "EMPTY";
  std::string covariance_file = "EMPTY";

  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"covariance", 	required_argument,0,'c'},
      {"number", 		required_argument,	0,'n'},
      {"signal", 		required_argument,	0,'s'},
      {"background", 	required_argument,	0,'b'},
      {"tag", 	    required_argument,	0,'t'},
      {"poisson", no_argument,0,'p'},
      {"help",no_argument,0,'h'},
      {0,			no_argument, 		0,  0},
    };

  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:n:s:b:c:t:ph", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 's':
	  signal_file = optarg;
	  break;
	case 'b':
	  background_file = optarg;
	  break;
	case 't':
	  tag = optarg;
	  break;
	case 'c':
      covariance_file = optarg;
	  break;
	case 'n':
	  num_MC_events = (int)strtod(optarg,NULL);
	  break;
    case 'p':
      sample_from_covariance = false;
      break;
	case '?':
	case 'h':
				std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"sbnfit_lee_frequentist_study allows for the simple hypothesis testing."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"--- Required arguments: ---"<<std::endl;
				std::cout<<"\t-x\t--xml\t\t\tInput configuration .xml file for SBNconfig"<<std::endl;
				std::cout<<"\t-t\t--tag\t\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"\t-b\t--background\t\tInput background only SBNspec.root file"<<std::endl;
	            std::cout<<"\t-c\t--covariance\t\tInput Fractional Covariance Matrix SBNcovar.root file"<<std::endl;
				std::cout<<"--- Optional arguments: ---"<<std::endl;
	            std::cout<<"\t-n\t--number\t\tNumber of MC events for frequentist studies (default 100k)"<<std::endl;
                std::cout<<"\t-p\t--poisson\t\tUse Stats-only covariance matrix with Poissonian statistics rather than use covariance matrix"<<std::endl;
				std::cout<<"\t-p\t--printall\t\tRuns in BONUS print mode, making individual spectra plots for ALLVariations. (warning can take a while!) "<<std::endl;
                std::cout<<"\t-h\t--help\t\t\tThis help menu."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
        return 0;	
	}
    }
  if(signal_file =="EMPTY"){
      std::cout<<"Error! You must enter a signal root file with the  `--signal  XX.SBNspec.root` or `-s XX.SBNspec.root` flags "<<std::endl;
      std::cout<<"Error! Run `--help` or `-h`  for more details."<<std::endl;
      return 1;
  }
   if(background_file =="EMPTY"){
      std::cout<<"Error! You must enter a background root file with the  `--background  XX.SBNspec.root` or `-b XX.SBNspec.root`  flags "<<std::endl;
      std::cout<<"Error! Run `--help` or `-h`  for more details."<<std::endl;
      return 1;
  }
  if(covariance_file =="EMPTY" && sample_from_covariance){
      std::cout<<"Error! You must enter a covariance root file with the  `--covariance  XX.SBNcovar.root` or `-c XX.SBNcovar.root`. "<<std::endl;
      std::cout<<"Error! Run `--help` or `-h`  for more details."<<std::endl;
      return 1;
  }
  


  std::cout<<"Loading signal file : "<<signal_file<<" with xml "<<xml<<std::endl;
  SBNspec sig(signal_file,xml);
  
  std::cout<<"Loading background file : "<<background_file<<" with xml "<<xml<<std::endl;
  SBNspec bkg(background_file,xml);

  std::cout<<"Loading fractional covariance matrix from "<<covariance_file<<std::endl;
  
  TFile * fsys = new TFile(covariance_file.c_str(),"read");
  TMatrixD * cov = (TMatrixD*)fsys->Get("frac_covariance");

  if(sample_from_covariance){
      SBNcls cls_factory(&bkg, &sig,*cov);
      cls_factory.SetSampleCovariance();
      cls_factory.CalcCLS(num_MC_events, tag);
  }else{
      SBNcls cls_factory(&bkg, &sig);
      cls_factory.CalcCLS(num_MC_events, tag);
  }

  return 0;
}
