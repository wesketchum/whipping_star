#include "SBNconditional.h"
using namespace sbn;

std::vector<TMatrixD> sbn::splitCovariance(TMatrixD & input, int start_cons_point){ 
        std::cout<<"SPLOT "<<input.GetNcols()<<std::endl;
        std::cout<<"sig sig "<<0<<" "<<start_cons_point-1<<" "<<0<<" "<<start_cons_point-1<<std::endl;
        TMatrixD sig_sig = input.GetSub(0,start_cons_point-1,0,start_cons_point-1); 
        TMatrixD cons_cons = input.GetSub(start_cons_point, input.GetNrows()-1 ,start_cons_point,input.GetNrows()-1); 
        TMatrixD sig_cons = input.GetSub(start_cons_point,input.GetNrows()-1,0,start_cons_point-1); 
        TMatrixD cons_sig = input.GetSub(0,start_cons_point-1,start_cons_point, input.GetNrows()-1); 


        std::vector<TMatrixD> ans = {sig_sig,sig_cons,cons_sig,cons_cons};
        return ans;
}


TMatrixD sbn::getConstrainedCovariance(std::vector<TMatrixD>& v_mat){
        
        TMatrixD cons_invert = v_mat.back();
        cons_invert.Invert();

        TMatrixD ans = v_mat.front()- v_mat[2]*cons_invert*v_mat[1];
        return ans;
}


