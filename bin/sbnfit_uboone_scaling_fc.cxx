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

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:r:p:sh", longopts, &index);

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
    mygrid.AddDimension("scale", 0.001, 8, 0.05);//0.1

    //Print the grid interesting bits
    mygrid.Print();
    SBNfeld myfeld(mygrid,tag,xml);


    if(mode_option == "feldman"){
    
        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        //myfeld.SetFractionalCovarianceMatrix(Msys);
        //myfeld.SetStatOnly();
        myfeld.m_subchannel_to_scale = "nu_uBooNE_1g1p_ncdelta";

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root","nu_uBooNE_1g1p_ncdelta",0.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        myfeld.FullFeldmanCousins();

    }else if(mode_option == "test"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        myfeld.LoadBackgroundSpectrum();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();

        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        myfeld.GlobalScan();

    }else if(mode_option == "scalescan"){


        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }
        myfeld.m_subchannel_to_scale = "nu_uBooNE_1g1p_ncdelta";

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root","nu_uBooNE_1g1p_ncdelta",0.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        TFile *f = new TFile("scan.root","recreate");
        f->cd();
        std::cout<<"Starting to peform a globalScan analysis"<<std::endl;
        std::vector<std::vector<double>> vec_grid = mygrid.GetGrid();
        TH2D * f_plot = new TH2D("f_plot","f_plot",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);
        TH2D * f_prob_plot = new TH2D("f_prob_plot","f_prob_plot",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);
        TH2D * f_prob_c_plot = new TH2D("f_prob_c_plot","f_prob_c_plot",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);
    
        std::cout<<"Makeing FC Maps!"<<std::endl;
        std::vector<TGraph*> tgr = myfeld.MakeFCMaps("SBNfeld_output_TEST.root");

        std::cout<<"MPrinting stuff"<<std::endl;
        for(int i=0; i< vec_grid.size(); i++){
            std::cout<<"On Grid Point "<<i<<" Scale_Factor: "<<vec_grid[i][0]<<std::endl;
            std::vector<double> vals = myfeld.GlobalScan2(i);
            for(int j=0; j < vec_grid.size(); j++){
                    f_plot->SetBinContent(j,i,vals[j]);
                    f_prob_plot->SetBinContent(j,i,tgr[i]->Eval(vals[j]));
                    f_prob_c_plot->SetBinContent(j,i,sqrt(2)*TMath::ErfInverse(1-tgr[i]->Eval(vals[j])));
                    std::cout<<"ScaleScan True: "<<vec_grid[i][0]<<" reco "<<vec_grid[j][0]<<" val "<<vals[j]<<" prob "<<tgr[i]->Eval(vals[j])<<std::endl;
            }
        }
        f->cd();
        f_plot->Write();
        f_prob_plot->Write();
        f_prob_c_plot->Write();
        TCanvas *c = new TCanvas("hope");
        c->cd();

        double contours1[1]; contours1[0]=1.0;
        double contours2[1]; contours2[0]=4.0;
        double contours3[1]; contours3[0]=9.0;
        double contours4[1]; contours4[0]=16.0; 

        f_plot->SetLineColor(kRed);
        f_plot->SetContour(1,contours1);
        f_plot->DrawCopy("cont3");

        f_plot->SetLineColor(kBlue);
        f_plot->SetContour(1,contours2);
        f_plot->DrawCopy("cont3 same");

        f_plot->SetLineColor(kGreen);
        f_plot->SetContour(1,contours3);
        f_plot->DrawCopy("cont3 same");

        f_plot->SetLineColor(kMagenta);
        f_plot->SetContour(1,contours4);
        f_plot->Draw("cont3 same");

        f_plot->GetXaxis()->SetTitle("Reconstructed Cross-Section (xSM)");
        f_plot->GetYaxis()->SetTitle("True Cross-Section (xSM)");


        c->Write();
        c->SaveAs("scan.png","png");
        f->Close();


    }
    else if(mode_option == "singlescan"){


        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }
        myfeld.m_subchannel_to_scale = "nu_uBooNE_1g1p_ncdelta";

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root","nu_uBooNE_1g1p_ncdelta",0.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout<<"Starting to peform a globalScan analysis"<<std::endl;
        std::vector<std::vector<double>> vec_grid = mygrid.GetGrid();

        {
            int i = grid_pt;

        std::cout<<"Makeing FC Maps!"<<std::endl;
        std::vector<TGraph*> tgr = myfeld.MakeFCMaps("SBNfeld_output_TEST.root",i);
            std::cout<<"On Grid Point "<<i<<" Scale_Factor: "<<vec_grid[i][0]<<std::endl;
            std::vector<double> vals = myfeld.GlobalScan2(i);
            for(int j=0; j < vec_grid.size(); j++){
                       std::cout<<"SimpleScan True: "<<vec_grid[i][0]<<" reco "<<vec_grid[j][0]<<" val "<<vals[j]<<" prob "<<tgr[i]->Eval(vals[j])<<std::endl;
                       //std::cout<<"SimpleScan: "<<j<<" "<<vec_grid[j][0]<<" "<<vals[j]<<std::endl;
            }
        }

    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

