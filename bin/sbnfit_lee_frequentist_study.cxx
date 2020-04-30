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
    bool sample_from_collapsed = false;
    bool remove_correlations = false;
    bool stats_only = false;
    int num_MC_events = 100000;
    bool use_cnp = false;
    int which_mode = 1;

    std::string tag = "TEST";

    std::string signal_file = "EMPTY";
    std::string background_file = "EMPTY";
    std::string covariance_file = "EMPTY";

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool zero_off_diag = false;
    bool tester=false;
    double epsilon = 1e-12;

    std::string legends;

    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"covariance", 	required_argument,0,'c'},
        {"collapse", 	required_argument,	0,'j'},
        {"signal", 		required_argument,	0,'s'},
        {"mode", 	    	required_argument,	0,'m'},
        {"background", 	required_argument,	0,'b'},
        {"tag", 	    required_argument,	0,'t'},
        {"epsilon", required_argument,0,'e'},
        {"legend",required_argument,0,'l'},
        {"cnp",no_argument,0,'a'},
        {"zero",no_argument,0,'z'},
        {"tester",no_argument,0,'k'},
        {"poisson", no_argument,0,'p'},
        {"flat", required_argument,0,'f'},
        {"help",no_argument,0,'h'},
        {0,			no_argument, 		0,  0},
    };

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "m:a:x:n:s:e:b:c:l:f:t:pjkzh", longopts, &index);

        switch(iarg)
        {
            case 'k':
                tester = true; 
                break;
            case 'z':
                remove_correlations = true; 
                break;
            case 'j':
                sample_from_collapsed = true;
                break;
            case 'm':
                which_mode = (int)strtod(optarg,NULL);
                break;
            case 'l':
                legends = optarg;
                break;
            case 'x':
                xml = optarg;
                break;
            case 's':
                signal_file = optarg;
                break;
            case 'b':
                background_file = optarg;
                break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
            case 'e':
                epsilon = (double)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'a':
                use_cnp = true;
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
                std::cout<<"\t-c\t--covariance\t\tInput Fractional Covariance Matrix SBNcovar.root file. If not passed, defaults to stats only!"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-j\t--collapse\t\tSample from collapsed rather than full covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t-f\t--flat\t\tAdd a flat percent systematic to fractional covariance matrix (all channels) (default false, pass in percent, i.e 5.0 for 5\% experimental)"<<std::endl;
                std::cout<<"\t-z\t--zero\t\tZero out all off diagonal elements of the systematics covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t-e\t--epsilon\t\tEpsilon tolerance by which to add back to diagonal of covariance matrix if determinant is 0 (default 1e-12)"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of MC events for frequentist studies (default 100k)"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tMode for test statistics 0: absolute chi^2, 1: delta chi^2 (default Delta Chi| obsolete, runs all concurrently)"<<std::endl;
                std::cout<<"\t-p\t--poisson\t\tUse Poissonian draws for pseudo experiments instead of from covariance matrix"<<std::endl;
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
    if(covariance_file =="EMPTY"){
        std::cout<<"Note! No covariance root file with the  `--covariance  XX.SBNcovar.root` or `-c XX.SBNcovar.root`. was passed, running in stats only mode. "<<std::endl;
        stats_only = true;
        sample_from_covariance = false;
    }


    std::cout<<"Loading signal file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec sig(signal_file,xml);

    std::cout<<"Loading background file : "<<background_file<<" with xml "<<xml<<std::endl;
    SBNspec bkg(background_file,xml);

    std::cout<<"Legends are being set to "<<legends.size()<<std::endl;

    std::cout<<"Loading fractional covariance matrix from "<<covariance_file<<std::endl;

    TFile * fsys;
    TMatrixD * cov;
    
    if(!stats_only){
        fsys = new TFile(covariance_file.c_str(),"read");
        cov = (TMatrixD*)fsys->Get("frac_covariance");
    }


    TMatrixD frac_flat_matrix(bkg.num_bins_total, bkg.num_bins_total);

    if(bool_flat_det_sys){
	    std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
	    frac_flat_matrix.ResizeTo(bkg.num_bins_total,bkg.num_bins_total);
	    frac_flat_matrix.Zero();//(bkg.num_bins_total,bkg.num_bins_total);
	    for(int i=0 ; i< bkg.num_bins_total; i++){
               frac_flat_matrix(i,i)=flat_det_sys_percent*flat_det_sys_percent/10000.;
        }
        std::cout<<"Just Before"<<std::endl;
        (*cov) = (*cov)+(frac_flat_matrix);
    }


    if(remove_correlations){
        std::cout<<"WARNING! We are running in   `Remove All Off Diagional Covariances/Correlations Mode` make sure this is what you want. "<<std::endl;
        for(int i=0; i<bkg.num_bins_total;i++){ 
            for(int j=0; j<bkg.num_bins_total;j++){ 
                if(i==j)continue;
                (*cov)(i,j) =0.0;
            }
        }
    }


    if(!stats_only){
        std::cout<<"Not running in stats only mode"<<std::endl;
        SBNcls cls_factory(&bkg, &sig,*cov);
        cls_factory.SetTolerance(epsilon);
        if(sample_from_collapsed)  cls_factory.SetSampleFromCollapsed();
        if(sample_from_covariance) cls_factory.SetSampleCovariance();
        cls_factory.setMode(which_mode);
        if(tester){cls_factory.runConstraintTest();return 0;}
        cls_factory.CalcCLS(num_MC_events, tag);
    }else{
        SBNcls cls_factory(&bkg, &sig);
        cls_factory.SetTolerance(epsilon);
        if(sample_from_collapsed)  cls_factory.SetSampleFromCollapsed();
        cls_factory.setMode(which_mode);
        if(tester){cls_factory.runConstraintTest();return 0;}
        cls_factory.CalcCLS(num_MC_events, tag);
    }

    return 0;
}
