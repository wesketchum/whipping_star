#include "SBNfit.h"
using namespace sbn;



SBNfit::SBNfit(SBNspec inBk, SBNspec inSg, TMatrixD mat, int npar) : SBNchi(inBk, mat), signal_osc_spectrum(inSg), num_params(npar) {

	for(int i =0; i< num_params; i++){
		f_is_fixed.push_back(0);
		f_param_names.push_back("");
		f_initial_values.push_back(0.5);
		f_upper_values.push_back(1);
		f_lower_values.push_back(0);
		f_step_sizes.push_back(0.01);
	}

	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
	f_minimizer_algo= "BFGS2";

	num_func_calls = 0;
}


SBNfit::SBNfit(SBNspec inBk, SBNspec inSg, int npar) : SBNchi(inBk), signal_osc_spectrum(inSg), num_params(npar) {


	for(int i =0; i< num_params; i++){
		f_is_fixed.push_back(0);
		f_param_names.push_back("");
		f_initial_values.push_back(0.5);
		f_upper_values.push_back(1);
		f_lower_values.push_back(0);
		f_step_sizes.push_back(0.01);
	}

	f_minimizer_mode ="GSLMultiMin"; //"GSLSimAn"
	f_minimizer_algo= "BFGS2";

	num_func_calls = 0;
}

int SBNfit::LoadSignal(SBNspec inSg){
	signal_osc_spectrum = inSg;
	return 0;
}

int SBNfit::InitializeNorm(std::vector< std::pair< std::string, int > > vecIn  ){
	vec_scales = vecIn;

	return 0;
}

double SBNfit::MinimizerCalcChi(const double * X){
	num_func_calls++;
	f_osc_spectrum = signal_osc_spectrum;

	for(auto& v: vec_scales){
		f_osc_spectrum.Scale(v.first, X[v.second] );
	}

	f_osc_spectrum.CollapseVector();

	last_calculated_chi =this->CalcChi(&f_osc_spectrum);

	return last_calculated_chi;

}

double SBNfit::Minimize(){
	num_func_calls=0;

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(f_minimizer_mode, f_minimizer_algo);

 	min->SetMaxIterations(20000);  // for GSL
	//min->SetTolerance(200000); //times 4 for normal
	min->SetPrintLevel(1);
	min->SetPrecision(0.001);//times 4 for normal


        ROOT::Math::Functor f( this, &SBNfit::MinimizerCalcChi, num_params);

   	min->SetFunction(f);

   	for(int i=0;i<num_params;i++){
		if(f_is_fixed[i]){
	   	min->SetFixedVariable(i,f_param_names[i],f_initial_values[i]);
	} else {
   		min->SetLimitedVariable(i,f_param_names[i],f_initial_values[i], f_step_sizes[i], f_lower_values[i],f_upper_values[i]);
	}

   	}

  	min->Minimize();

  	const double *xs = min->X();

    bf_chi= MinimizerCalcChi(xs);;
    bf_params = xs;
    return bf_chi;

}

/****************************************************
 ***		Some initial Setup things
 * *************************************************/

int SBNfit::SetMethod(std::string mode, std::string algo){
	f_minimizer_mode =mode; //"GSLSimAn"
	f_minimizer_algo= algo;


return 0;
}


int SBNfit::SetInitialValues(std::vector<double> inv){
	for(int i = 0; i< num_params; i++){
		f_initial_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::SetInitialValues(double in){
	for(int i = 0; i< num_params; i++){
		f_initial_values[i]=in;
	}
	return 0;
}



int SBNfit::SetUpperValues(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_upper_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::SetUpperValues(double in){
	for(int i = 0; i< num_params; i++){
		f_upper_values[i]=in;
	}
	return 0;
}



int SBNfit::SetLowerValues(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_lower_values[i]=inv[i];
	}
	return 0;
}
int SBNfit::SetLowerValues(double in){
	for(int i = 0; i< num_params; i++){
		f_lower_values[i]=in;
	}
	return 0;
}

int SBNfit::SetStepSizes(std::vector<double>  inv){
	for(int i = 0; i< num_params; i++){
		f_step_sizes[i]=inv[i];
	}
	return 0;
}
int SBNfit::SetStepSizes(double in){
	for(int i = 0; i< num_params; i++){
		f_step_sizes[i]=in;
	}
	return 0;
}



int SBNfit::SetFixed(std::vector<int>  inv){
	for(int i = 0; i< num_params; i++){
		f_is_fixed[i]=inv[i];
	}
	return 0;
}
int SBNfit::SetFixed(int in){
	for(int i = 0; i< num_params; i++){
		f_is_fixed[i]=in;
	}
	return 0;
}


int SBNfit::SetNames(std::vector<std::string> inv){
	for(int i = 0; i< num_params; i++){
		f_param_names[i]=inv[i];
	}
	return 0;
}
