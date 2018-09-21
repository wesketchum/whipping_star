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

#include "SBNsvd.h"

#include <stdio.h>
#include <cuda.h>
#include "magma_v2.h"
#include "magma_lapack.h"
//#include "magma_dutil.cpp"


#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main ( int argc , char ** argv ){
    
    magma_init();
    int size[10] = {50,100,200,400,800,1600,3200,6400,12500,25000};


    for(int i=0; i< 10; i++){
    magma_queue_t queue = NULL ;
    magma_int_t dev =0;
    magma_queue_create (dev ,& queue );
    double gpu_time , * dwork ; // dwork - workspace
    magma_int_t ldwork ; // size of dwork
    magma_int_t *piv , info ; // piv - array of indices of inter -
    magma_int_t m = size[i]; // changed rows ; a - mxm matrix
    magma_int_t mm=m*m; // size of a, r, c
    double *a; // a- mxm matrix on the host


    magma_int_t ione = 1;
    magma_int_t ISEED [4] = {0 ,0 ,0 ,1 }; // seed
    magma_int_t err;
    const double alpha = 1.0; // alpha =1
    const double beta = 0.0; // beta =0
    ldwork = m * magma_get_dgetri_nb ( m ); 

    {
    double *d_a ; // d_a - mxm matrix a on the device
    double *d_r ; // d_r - mxm matrix r on the device
    double *d_c ; // d_c - mxm matrix c on the device
    // allocate matrices
    err = magma_dmalloc_cpu (&a , mm ); // host memory for a
    err = magma_dmalloc (&d_a , mm ); // device memory for a
    err = magma_dmalloc (&d_r , mm ); // device memory for r
    err = magma_dmalloc (&d_c , mm ); // devicie memory for c

    err = magma_dmalloc ( &dwork , ldwork ); // dev. mem. for ldwork
    piv =( magma_int_t *) malloc (m* sizeof ( magma_int_t )); // host mem.
     // generate random matrix a // for piv
     //
    lapackf77_dlarnv (& ione ,ISEED ,&mm ,a); // randomize a
    magma_dsetmatrix ( m, m, a,m, d_a ,m, queue ); // copy a -> d_a
    magmablas_dlacpy ( MagmaFull ,m,m,d_a ,m,d_r ,m, queue ); // d_a - >d_r

    // find the inverse matrix : d_a *X=I using the LU factorization
    // // with partial pivoting and row interchanges computed by
    // magma_dgetrf_gpu ; row i is interchanged with row piv (i);
    // d_a -mxm matrix ; d_a is overwritten by the inverse
     
     gpu_time = magma_sync_wtime ( NULL );

     //std::cout<<"About to run dgetrf"<<std::endl;
     magma_dgetrf_gpu( m, m, d_a, m, piv, &info);
     magma_dgetri_gpu(m,d_a,m,piv,dwork,ldwork,&info);
     //std::cout<<"Ran to run dgetrf"<<std::endl;
     gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
     magma_dgemm ( MagmaNoTrans , MagmaNoTrans ,m,m,m,alpha ,d_a ,m,  d_r ,m,beta ,d_c ,m, queue ); // multiply a^ -1*a
     printf (" magma_dgetrf_gpu + magma_dgetri_gpu time :%d %7.5f sec  \n",size[i],gpu_time );
     magma_dgetmatrix ( m, m, d_c , m, a, m, queue ); // copy d_c - >a
     //printf (" upper left corner of a^ -1*a:\n");
    //magma_dprint ( 4 , 4 , a, m ); // part of a^ -1*a
     
     free (a); // free host memory
     magma_free (d_a ); // free device memory
     magma_free (d_r ); // free device memory
     magma_free (d_c ); // free device memory
     magma_queue_destroy ( queue ); // destroy queue

     


    }
     free ( piv ); // free host memory

    }


     magma_finalize (); // finalize Magma
     return 0;
     }
     // magma_dgetrf_gpu + magma_dgetri_gpu time : 4.79694 sec .
     // upper left corner of a^ -1*a:
     //[
     // 1.0000 -0.0000 -0.0000 0.0000
     // 0.0000 1.0000 -0.0000 -0.0000
     // -0.0000 0.0000 1.0000 0.0000
     // -0.0000 0.0000 -0.0000 1.0000
     // ];
