#ifndef SBNFIT_H_
#define SBNFIT_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

#include "SBNconfig.h"
#include "SBNspec.h"
#include "SBNchi.h"

#include <TH1D.h>
#include <string>
#include <TF1.h>
#include <TMatrixT.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
namespace sbn{


class SBNfit : public SBNchi {
	
	protected:

	std::vector<double> f_initial_values;
	std::vector<double> f_upper_values;
	std::vector<double> f_lower_values;
	std::vector<double> f_step_sizes;
	std::vector<int>    f_is_fixed;
	std::vector<std::string> f_param_names;

	std::string f_minimizer_mode;
	std::string f_minimizer_algo;


	SBNspec f_osc_spectrum;

	public:
	//using SBNchi::SBNchi;
	SBNfit(SBNspec bk, SBNspec sk,int npa);
	
	SBNfit(SBNspec inBk, SBNspec inSg, TMatrixD mat, int npar);



	SBNspec signal_osc_spectrum;

	int num_params;
	std::vector<std::pair<std::string, int>> vec_scales;
	double * fX;

	double bf_chi;
	const double * bf_params;
	
	int num_func_calls;

	

	int LoadSignal(SBNspec);

	int InitializeNorm(std::vector< std::pair<std::string, int>> );
	virtual double MinimizerCalcChi(const double * X);
	double Minimize();	

	//ROOT::Math::Minimizer* min ;     

	int SetMethod(std::string, std::string);

	int SetInitialValues(std::vector<double>);
	int SetInitialValues(double);
	
	int SetUpperValues(std::vector<double>);
	int SetUpperValues(double);

	int SetLowerValues(std::vector<double>);
	int SetLowerValues(double);
	
	int SetStepSizes(std::vector<double>);
	int SetStepSizes(double);
	
	int SetFixed(std::vector<int>);
	int SetFixed(int);
	
	int SetNames(std::vector<std::string>);

};


};
#endif
