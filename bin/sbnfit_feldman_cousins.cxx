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
        {"flat",        required_argument, 0 ,'f'},
        {"randomseed",        required_argument, 0 ,'r'},
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
    int number = 2500;
    int grid_pt = 0;
    double random_number_seed = -1;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:r:p:f:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'n':
                number = (int)strtod(optarg,NULL);
                break;
            case 'p':
                grid_pt = (int)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
            case 'm':
                mode_option = optarg;
                break;
            case 'r':
                random_number_seed = (double)strtod(optarg,NULL);
                std::cout<<"Reading in random seed argument: "<<random_number_seed<<std::endl;
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
                std::cout<<"\t-m\t--mode\t\tWhat mode you want to run in. Arguments are:"<<std::endl;
                std::cout<<"\t\t\t--\t gen : Generate the preoscillated spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t\t\t--\t genbkg : Generate a background only spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t\t\t--\t feldman : Perform a full feldman cousins analysis"<<std::endl;  
                std::cout<<"\t-p\t--point\t\tWhat Grid Point to run over. -1 for a Full run (WARNING takes forever) [Defaults to 0]"<<std::endl;
                std::cout<<"\t\t\t--\t test : Just a testbed. Can be ignored"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of pseudo-experiments to simulate (default 2500)"<<std::endl; 
                std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
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

    if(tag == "NuMuDis"){	
       //grid for numu disappearance
       mygrid.AddDimension("m4", -1, 1.05, 0.05);//0.05
       mygrid.AddFixedDimension("ue4", 0);
       mygrid.AddDimension("um4",-1.2, -0.5, 0.01); //for NuMuAllowed
       //mygrid.AddDimension("um4",-2.0, 0.05, 0.05); //0.05
    }else{
      //grid for nue appearance
      mygrid.AddDimension("m4", -1.0, 1.05, 0.05);   //0.1
      mygrid.AddDimension("ue4", -2.3, 0.05, 0.05); //0.1
      mygrid.AddFixedDimension("um4",0.0);         //0.05
    }


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

    }else if(mode_option == "feldman"){

        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        std::cout <<"DONE loading precomputed spectra at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
        myfeld.LoadBackgroundSpectrum();

        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        if(grid_pt ==-1){
            std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
            myfeld.FullFeldmanCousins();
        }else if(grid_pt>=0){
             std::cout<<"Beginning to peform Single Grid PointFeldmanCousins analysis on pt: "<<grid_pt<<std::endl;
             myfeld.PointFeldmanCousins((size_t)grid_pt);
        }

    }else if(mode_option == "test"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
 
          if(number >= 0) {
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance_"+std::to_string(number));
          }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
          }
        }

        if(bool_flat_det_sys){
            myfeld.AddFlatDetSystematic(flat_det_sys_percent);
        }


        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        myfeld.LoadBackgroundSpectrum();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();

        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        //myfeld.GlobalScan(884);
        myfeld.GlobalScan(1503);

    }else if(mode_option == "plot"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        myfeld.SetEmptyFractionalCovarianceMatrix();
        myfeld.SetStatOnly();
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.LoadPreOscillatedSpectrum(1503);

    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

