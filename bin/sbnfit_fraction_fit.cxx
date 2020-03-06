#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TList.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TROOT.h"

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
        {"montecarlo", 		required_argument,	0, 'm'},
        {"data", 		required_argument,	0, 'd'},
	{"interpolation",       required_argument,      0, 'i'},
	{"covarmatrix",         required_argument,      0, 'c'},
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
    bool bool_stat_only = false;
    int interpolation_number = -99;  //number of points for chi2 value interpolation
    double random_number_seed = -1;

    bool input_data = false;
    std::string data_filename;
    std::string mc_filename;
    std::string covmatrix_file;  //root file containing covariance matrix

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "d:x:m:r:p:i:f:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
                break;
	    case 'm':
		mc_filename = optarg;
		break;
	    case 'i':
		interpolation_number = (int)strtod(optarg, NULL);
		break;
	    case 'c':
		covmatrix_file = optarg;
		break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
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
                std::cout<<"sbnfit_fraction_fit is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-d\t--data\t\tInput observed data for a global scan"<<std::endl;
                std::cout<<"\t-m\t--montecarlo\t\tInput monte carlo for global scan"<<std::endl;
		std::cout<<"\t-i\t--interpolation\t\tInput number of points for interpolation"<< std::endl;
		std::cout<<"\t-c\t--covariance matrix\t\tInput syst covariance matrix for global scan"<< std::endl;
                std::cout<<"\t-f\t--flat\t\t Input flat systematic fractional covariance matrix"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
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

    std::cout<<"Begining Fraction Fit for tag: "<<tag<<std::endl;

    NGrid mygrid;

    //now only available for 2 subchannels only
    mygrid.AddDimension("BNBoverlay", -0.55, 0.55, 0.05);   //0.1 FULL
    mygrid.AddDimension("ncpi0overlay", -1.0, 1.05, 0.05);   //0.1 FULL


    std::cout << "Fraction Fit|| "<<tag << "\tStart initializing MC and data spectrum" << std::endl;

    //initialize the MC spectrum
    SBNspec mc_spec(mc_filename, xml);

    //initlaize the data spectrum
    SBNspec data_spec(data_filename, xml);
    data_spec.CollapseVector();  //collapse full vector



    std::cout << "Fraction Fit||" << tag << "\tInitialize fractional systematric covariance matrix" << std::endl;
    //initialize covariance matrix
    TMatrixT<double> frac_syst_matrix(mc_spec.num_bins_total, mc_spec.num_bins_total);


    if(bool_stat_only){
    	frac_syst_matrix.Zero();
        std::cout<<"RUNNING Stat Only!"<<std::endl;
    }else if(bool_flat_det_sys){
	std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
	frac_syst_matrix.Zero();
	//set up flat fractional syst covariance
	for(int i=0 ; i< mc_spec.num_bins_total; i++)
		for(int j=0;j<mc_spec.num_bins_total; j++)
			frac_syst_matrix(i,j)=flat_det_sys_percent*flat_det_sys_percent/10000.;
    }
    else{
	std::cout<< "RUNNING with Systematics!" << std::endl;
	frac_syst_matrix.Zero();

	std::cout << "Open covariance matrix root file: " << covmatrix_file << std::endl;
	TFile* f_covar = new TFile(covmatrix_file.c_str(), "read");
	TMatrixT<double>* m_temp;
	m_temp = (TMatrixD*)f_covar->Get("frac_matrix_name");
	frac_syst_matrix = *m_temp;
	f_covar->Close();	
    }



   std::cout<< "Fraction Fit||" << tag<< "\tGrab info of the grid"<< std::endl;
   //check grid size
   std::vector<std::vector<double>> grid = mygrid.GetGrid();
   if(grid.size() != mygrid.f_num_total_points){
	std::cout <<  "the number of points don't match: something wrong with the grid setup!!" << std::endl;
	return 1;
   }

   //collect the name of dimensions: subchannels you want to vary; and the range
   std::vector<std::string> dimension_name;
   double range_x_low = mygrid.f_dimensions.at(0).f_min;
   double range_x_up = mygrid.f_dimensions.at(0).f_max;
   double range_y_low = mygrid.f_dimensions.at(1).f_min;
   double range_y_up = mygrid.f_dimensions.at(1).f_max;
   int nbin_x = mygrid.f_dimensions.at(0).f_N;  //number of point in x axis
   int nbin_y = mygrid.f_dimensions.at(1).f_N;  
 
   dimension_name.clear(); 
   for(int i=0; i< mygrid.f_num_dimensions ; i++){
	dimension_name.push_back(mygrid.f_dimensions.at(i).f_name);
   }



   TFile* f_output = new TFile("chi_contour.root", "recreate");
   TH2D* h_chi2_raw = new TH2D("h_chi2_raw", "h_chi2_raw", nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   TH2D* h_chi2_delta = new TH2D("h_chi2_delta", "h_chi2_delta", nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   TH2D* h_chi2_inter = new TH2D("h_chi2_interpolation", "h_chi2_interpolation", nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   
   //*********************loop over grid points, calc chi square********************************


   std::cout << "Fraction Fit||"<< tag<< "\tStart GLOBAL SCAN" <<std::endl;
   std::vector<double> chi;  //vector to save chi square values.
   chi.reserve(grid.size());  //reserve the memory

   for(int i =0; i< grid.size() ;i++){

	//set a temperary SBNspec, assume there is already a MC CV root file with corresponding sequence defined in xml	
	SBNspec spec_temp(mc_filename, xml, false);
	
	//access grid point
	std::vector<double> point = grid[i];

	//check dimension point
	if(point.size() != dimension_name.size()){
		std::cout << "dimension doesn't match: something wrong with grid setup!" << std::endl;
		return 1;
	}

	//scale chosen subchannels
	for(int j=0; j< point.size(); j++ ){
		spec_temp.Scale(dimension_name[j], point[j]);
	}

	std::cout << "Point " << i << " FULL VECTOR: ";
	for(int j =0; j<spec_temp.full_vector.size(); j++) std::cout << " " << spec_temp.full_vector.at(j) ;
	std::cout << std::endl;	
	
	std::cout << "Point " << i << " COLLAPSED VECTOR: ";
	for(int j =0; j<spec_temp.collapsed_vector.size(); j++) std::cout << " " << spec_temp.collapsed_vector.at(j) ;
	std::cout << std::endl;	

	//setup SBNchi with spectrum expected at this grid point
	SBNchi chi_temp(spec_temp, frac_syst_matrix);
	//calculate chi2
	chi.push_back(chi_temp.CalcChi(data_spec.collapsed_vector));
	h_chi2_raw->Fill(point[0], point[1], chi.back());

	//print out chi square value
	std::cout << "CHI2: " << i << " " << chi.back();
	for(int j=0; j<point.size() ; j++) std::cout << " " << point[j];
	std::cout << std::endl;

   }
   
   f_output->cd();
   h_chi2_raw->Write();   

   double chi_min=DBL_MAX; // minimum of chi2.
   chi_min = *std::min_element(chi.begin(), chi.end());

   std::cout << "Print out delta chi square for grid points now" << std::endl;
   //delta chi2;
   for(int i=0; i< grid.size(); i++){
	 chi[i] -= chi_min;
	 std::vector<double> point = grid[i];
	 h_chi2_delta->Fill(point[0], point[1], chi[i]);	 
	 std::cout << "DeltaChi: " << i << " " << chi[i];
	 for(int j =0;j <point.size(); j++) std::cout << " " << point[j];
	 std::cout << std::endl;
   }
   f_output->cd();
   h_chi2_delta->Write();



   
   //****************** START INTERPOLATION*********************************

   if(interpolation_number != -99){
	std::cout << "Fraction Fit||" << tag << "\tStart interpolation with number "<< interpolation_number <<std::endl;
   	const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	std::vector<double> grid_x;
	std::vector<double> grid_y;
	grid_x = mygrid.f_dimensions.at(0).f_points;  //initialize grid point
	grid_y = mygrid.f_dimensions.at(1).f_points;   //initialize grid point
	double x_step = ceil(fabs(range_x_low - range_x_up)/interpolation_number);
	double y_step = ceil(fabs(range_y_low - range_y_up)/interpolation_number);
  	gsl_spline2d *spline = gsl_spline2d_alloc(T, nbin_x, nbin_y);
  	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  	gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  	size_t i, j;

	/* initialize interpolation */
	gsl_spline2d_init(spline,  &grid_y[0], &grid_x[0], &chi[0] , nbin_y, nbin_x);

	/* interpolate N values in x and y and print out grid for plotting */
	for (i = 0; i < interpolation_number; ++i){
	      double xi = range_x_low + x_step*i;;
	      for (j = 0; j < interpolation_number; ++j)
		{
		  double yj = range_y_low + y_step*j;;
		  double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);
		  h_chi2_inter->Fill(xi, yj, zij);
		  std::cout << "Interpolated value: " << zij << " " << xi << " " << yj << std::endl; 
		}
	}

	  f_output->cd();
	  h_chi2_inter->Write();

	  gsl_spline2d_free(spline);
	  gsl_interp_accel_free(xacc);
	  gsl_interp_accel_free(yacc);
  }
  //*****************END OF INTERPOLATION*************************************************






  //****************draw contour, get range for fraction fit****************************
   std::vector<double> contour{2.30, 4.61,6.18,11.83}; // chi2 value for 2 dof with 1 sigma, 90%, 2 sigma and 3 sigma confidence level
   TCanvas* c = new TCanvas("c", "c");
   h_chi2_inter->SetContour((int)contour.size(), &contour[0]);
   h_chi2_inter->Draw("CONT Z LIST");//"LIST" generates a list of TGraph for each contour
   c->Update();
   f_output->cd();
   c->Write("contour");

   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

   TH2D* hr = (TH2D*)h_chi2_delta->Clone();
   hr->Scale(0); // keep onlythe x and y axis

   TCanvas *v_canvas[conts->GetSize()];
   TGraph* c_graph=NULL;
   TGraph* copy_graph = NULL;
   TList* con_list=NULL;
   for(int i=0; i<conts->GetSize(); i++){
	con_list = (TList*)conts->At(i);   //a list of TGraph for i'th contour.
	std::cout << "Contour Z=" << contour[i] << " has more than 1 enclosed region! " << std::endl;
	std::cout << "Size of the list is " << con_list->GetSize() << std::endl;

	double x_min = DBL_MAX;
	double x_max = DBL_MIN;
	double y_min = DBL_MAX;
	double y_max= DBL_MIN;
	v_canvas[i] = new TCanvas(Form("c%d", i), Form("c%d", i));
        hr->Draw();
	c_graph = (TGraph*)con_list->First();  //grab the TGraph

	for(int j=0; j< con_list->GetSize() ; j++){
		copy_graph= (TGraph*)c_graph->Clone();
                copy_graph->Draw("C");
		std::cout << j << ": number of point " << copy_graph->GetN() << std::endl;
			
		// grab x ,y coordinate
		double x,y;
		for(int k =0; k< copy_graph->GetN(); k++){
			copy_graph->GetPoint(k, x, y);
			if(x < x_min) x_min = x;
			else if(x > x_max) x_max = x;

			if(y < y_min) y_min =y;
			else if(y > y_max) y_max = y;
		}
	
		//move to next graph
		c_graph=(TGraph*)con_list->After(c_graph);
	}

	v_canvas[i]->Update();
	f_output->cd();
	v_canvas[i]->Write(Form("contour_%f",contour[i]));

	std::cout << "Contour " << contour[i] << std::endl;
	std::cout << "range for x :" << x_min << "~" << x_max << std::endl;
	std::cout << "range for y : " << y_min << "~" << y_max << std::endl;
   }
  
   


    std::cout << "Actually done with everything! " << std::endl;
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

