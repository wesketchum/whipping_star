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

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
  std::string xml = "example.xml";
  int iarg = 0;
  opterr=1;
  int index;
  bool sample_from_covariance = false;
  int num_MC_events = 40000;

  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"covariance", 		no_argument,0,'c'},
      {"number", 		required_argument,	0,'n'},
      {0,			no_argument, 		0,  0},
    };

  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:n:c", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 'c':
	  sample_from_covariance = true;
	  break;
	case 'n':
	  num_MC_events = (int)strtod(optarg,NULL);
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  std::cout<<"\t-c\t--covariance\t\tSample from covariance matrix instead of Poisson"<<std::endl;
	  std::cout<<"\t-n\t--number\t\tNumber of MC events for frequentist studies (default 40k)"<<std::endl;
	  return 0;
	}
    }

  std::string tag = "EXAMPLEI";

    SBNspec sig("EXAMPLE1.SBNspec.root",xml);
    sig.Scale("leesignal",1.5);


    magma_init();
    int size[14] = {3,6,12,25,50,100,200,400,800,1600,3200,6400,12500,25000};

    for(int i=0; i< 12;i++){

    magma_queue_t queue = NULL ;
    magma_int_t dev =0;
    magma_queue_create (dev ,& queue );
    double gpu_time, gpu_time_mag, * dwork ; // dwork - workspace
    magma_int_t ldwork ; // size of dwork
    magma_int_t *piv , info ; // piv - array of indices of inter -
    magma_int_t m = size[i]; // changed rows ; a - mxm matrix
    int n_t = m;
    magma_int_t mm=m*m; // size of a, r, c
    double *a; // a- mxm matrix on the host


    magma_int_t ione = 1;
    magma_int_t ISEED [4] = {0 ,0 ,0 ,1 }; // seed
    magma_int_t err;
    const double alpha = 1.0; // alpha =1
    const double beta = 0.0; // beta =0
    ldwork = m * magma_get_dgetri_nb ( m ); 

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
     


    double** a_matrix = new double*[n_t];
    double** a_matrix2 = new double*[n_t];
    double** a_inverse = new double*[n_t];
    double** a_mult = new double*[n_t];

    for(int i=0; i < n_t; i++){
        a_matrix[i] = new double[n_t];
        a_matrix2[i] = new double[n_t];
        a_inverse[i] = new double[n_t];
        a_mult[i] = new double[n_t];
    }

    for(int i=0; i< n_t; i++){
    
        for(int j=0; j< n_t; j++){
            a_matrix[i][j] =  a[i*n_t+j];
            a_matrix2[i][j] = a[i*n_t+j];
            a_mult[i][j]=0.0;

            if(i==j){
                 a_inverse[i][j]= 1.0;
            }else{
                a_inverse[i][j] = 0.0;
                
            }
        }
    }

     gpu_time = magma_sync_wtime ( NULL );
     MatrixInverse(a_matrix, a_inverse, n_t);
     gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
    

     MatrixMult(a_matrix2, a_inverse, a_mult, n_t);


     gpu_time_mag = magma_sync_wtime ( NULL );
     //std::cout<<"About to run dgetrf"<<std::endl;
     magma_dgetrf_gpu( m, m, d_a, m, piv, &info);
     magma_dgetri_gpu(m,d_a,m,piv,dwork,ldwork,&info);
     //std::cout<<"Ran to run dgetrf"<<std::endl;
     gpu_time_mag = magma_sync_wtime ( NULL ) - gpu_time_mag ;
     magma_dgemm ( MagmaNoTrans , MagmaNoTrans ,m,m,m,alpha ,d_a ,m,  d_r ,m,beta ,d_c ,m, queue ); // multiply a^ -1*a
     magma_dgetmatrix ( m, m, d_c , m, a, m, queue ); // copy d_c - >a
     //printf (" upper left corner of a^ -1*a:\n");
    //magma_dprint ( 4 , 4 , a, m ); // part of a^ -1*a

     printf ("MINE :%d %7.5f sec : MAG :  %7.5f \n",size[i],gpu_time,gpu_time_mag );
     free (a); // free host memory
     magma_free (d_a ); // free device memory
     magma_free (d_r ); // free device memory
     magma_free (d_c ); // free device memory
     magma_queue_destroy ( queue ); // destroy queue
     free ( piv ); // free host memory


   for(int i=0; i < n_t; i++){
        delete[] a_matrix[i];  
        delete[] a_inverse[i];  
        delete[] a_mult[i];  
    }

    delete[] a_matrix;
    delete[] a_inverse;
    delete[] a_mult;
    

    }


     magma_finalize (); // finalize Magma
     

















    /*std::cout<<"\n"<<"InverOut\n";
   for(int i=0; i < n_t; i++){
   for(int j=0; j < n_t; j++){
       std::cout<<a_mult[i][j]<<"\t"; 
   }
   std::cout<<"\n";
   }
  */

  return 0;
}
