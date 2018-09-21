#ifndef SBNSVD_H_
#define SBNSVD_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNconfig.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"

#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


namespace sbn{





int MatrixInverse(double ** matrix, double ** inverse,int);
int MatrixMult(double ** a, double ** b, double **c, int n );

int SingleValueDecomposition(double ** matrix, double ** U, double**V, double *single_values );
int Pythag(double a, double b); 

int MatrixInverse(double ** matrix, double ** inverse, int n ){


    // #pragma acc parallel loop copy(matrix[:n][:n],inverse[:n][:n])
    for(int i=0; i<n; i++){

    int temp=i;


    for(int j=i+1; j<n; j++){
        if(matrix[j][i] > matrix[temp][i])    temp=j;
    }
      
    if(fabs(matrix[temp][i])<1e-14){
          //  std::cout<<"Aghr: "<<matrix[temp][i]<<std::endl;
    //        exit(0);
     }


        if(temp!=i){
              for(int k=0; k<n; k++){

                double tmp_matrix_k = matrix[i][k];
                matrix[i][k] = matrix[temp][k];
                matrix[temp][k] = tmp_matrix_k;

                double tmp_inverse_k = inverse[i][k];
                inverse[i][k] = inverse[temp][k];
                inverse[temp][k] = tmp_inverse_k;
            }
        }

        for(int j=0; j<n; j++){

            if(j!=i){
                double tmp = matrix[j][i];

                //k is rows again
                for(int k=0; k<n; k++){
                    matrix[j][k]  -=  tmp*(matrix[i][k]/matrix[i][i]);
                    inverse[j][k]  -=  tmp*(inverse[i][k]/matrix[i][i]);
                }

            }else{
                double tmp = matrix[j][i];

                for(int k=0; k<n; k++){
                    matrix[j][k] /=  tmp;
                    inverse[j][k] /=  tmp;
                }
            }
        }
    }
    return 0;
};



int MatrixMult(double ** a, double ** b, double **c, int n ){
    //v = a*b

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            c[i][j]=0;
            for(int k=0; k<n; k++){
                c[i][j] += a[i][k]*b[k][j];
            }

        }
    }

    return 0;
};

}
#endif
