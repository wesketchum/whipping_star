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
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_make_covariance.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "oscillate_example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"printall", 		no_argument, 		0, 'p'},
        {"stat", 		no_argument, 		0, 's'},
         {"number", 		required_argument,	0,'n'},
   {"tag", 		required_argument,	0, 't'},
        {"mode",        required_argument, 0 ,'m'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    std::string mode_option;
    bool bool_stat_only = false;
    int number = 0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:psh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
	case 'n':
	  number = (int)strtod(optarg,NULL);
	  break;

            case 't':
                tag = optarg;
                break;
            case 'm':
                mode_option = optarg;
                break;
            case 's':
                bool_stat_only = true;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_feldman_cousins is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
    }

    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"Begining FeldmanCousins for tag: "<<tag<<std::endl;

    NGrid mygrid;
    
 //   mygrid.AddDimension("m4", -1, 1.05, 0.05);//0.05
 //   mygrid.AddFixedDimension("ue4", 0);
 //   mygrid.AddDimension("um4",-2.0, 0.1, 0.1); //0.05

    //1 -1 -0.7
    //2 -0.7 -0.4
    //3 -0.4 -0.1
    //4 -0.1 0.2
    //5 0.2 0.5
    //6 0.5 0.8
    //7 0.8 1.1
    
    std::vector<double> low = {-1.0,-0.7,-0.4,-0.1,0.2,0.5,0.8};
    std::vector<double> hi = {-0.7,-0.4,-0.1,0.2,0.5,0.8,1.1};

    mygrid.AddDimension("m4", low[0], hi[6], 0.1);//0.05
    mygrid.AddDimension("ue4", -2.3, 0.1, 0.1);
    mygrid.AddFixedDimension("um4",0.0); //0.05



    //Print the grid interesting bits
    mygrid.Print();
    SBNfeld myfeld(mygrid,tag,xml);

    if(mode_option == "gen"){
        myfeld.GenerateOscillatedSpectra();
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.GenerateBackgroundSpectrum();

    }else if(mode_option == "genbkg"){
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.GenerateBackgroundSpectrum();

    }else if(mode_option == "fit"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        myfeld.LoadBackgroundSpectrum();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();

        std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        myfeld.FullFeldmanCousins();

    }else if(mode_option == "global"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        
        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }

        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        myfeld.LoadBackgroundSpectrum();


        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();

        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        myfeld.GlobalScan();

    }else if(mode_option == "test"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");




    }
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

