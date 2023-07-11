#ifndef BRANCHVARIABLE_H_
#define BRANCHVARIABLE_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <time.h>


struct BranchVariable{
	std::string name;
	std::string type;
	std::string associated_hist;

	bool oscillate;
	std::string true_param_name;
	std::string true_L_name;

	//int value_i;
	//int true_value_i;
	//int true_L_i;

	float value_f;
	float true_value_f;
	float true_L_f;

	double value_d;
	double true_value_d;
	double true_L_d;

	BranchVariable(std::string n, std::string t, std::string a) : name(n), type(t), associated_hist(a) {oscillate=false;};
	virtual void* GetValue(){};
	virtual void* GetTrueValue(){};
	virtual void* GetTrueL(){};

	int SetOscillate(bool inbool){ oscillate = inbool;};
	bool GetOscillate(){ return oscillate;};

};

/*struct BranchVariable_i: public BranchVariable{
	BranchVariable_i(std::string n, std::string t, std::string a) : BranchVariable(n,t,a) {value_i=0;true_value_i=0; true_L_i =0;};
	void* GetValue(){ return &value_i;}
	void* GetTrueValue(){ return &true_value_i;}
	void* GetTrueL(){ return &true_L_i;}
};
*/

struct BranchVariable_d: public BranchVariable{
	BranchVariable_d(std::string n, std::string t, std::string a) : BranchVariable(n,t,a) {value_d=0;true_value_d=0; true_L_d = 0;};
	void* GetValue(){ return &value_d;}
	void* GetTrueValue(){ return &true_value_d;}
	void* GetTrueL(){ return &true_L_d;}
};

struct BranchVariable_f: public BranchVariable{
	BranchVariable_f(std::string n, std::string t, std::string a) : BranchVariable(n,t,a) {value_f=0;true_value_f=0; true_L_f = 0;};
	void* GetValue(){ return &value_f;}
	void* GetTrueValue(){ return &true_value_f;}
	void* GetTrueL(){ return &true_L_f;}
};

#endif
