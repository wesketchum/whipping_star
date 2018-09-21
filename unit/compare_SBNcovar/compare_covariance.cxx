#include <iostream>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TMatrixD.h"
#include "TStopwatch.h"

TStopwatch global_watch, local_watch;

void StartWatch() {
   local_watch.Reset(); 
   local_watch.Start();
}
void StopWatch() { 
  local_watch.Stop();
  std::cout << "(" << local_watch.CpuTime() 
	    << "," << local_watch.RealTime() << ")" << std::endl;
}

void Compare(const TMatrixD* m1, const TMatrixD* m2) {
  static double EPS=1e-6;

  int nrows = m1->GetNrows();
  int ncols = m1->GetNcols();

  StartWatch();
  for(int i=0; i<nrows; ++i) {
    for(int j=0; j<ncols; ++j) {
      if (std::abs((*m1)(i,j) - (*m2)(i,j)) > EPS) {
	std::stringstream ss; 
	ss << "Difference" << std::endl;
	ss << "@(i,j)=("<<i<<","<<j<<") m1="<<(*m1)(i,j)<<" m2=" << (*m2)(i,j) << std::endl;
	std::exit(1);
	return;

      }
    }
  } 
  StopWatch();

}

int main(int argc, char** argv) {
  global_watch.Reset();
  global_watch.Start();

  if (argc != 3) {
    std::cout << std::endl;
    std::cout << "argv[1] = EXAMPLE1.SBNcovar.root" << std::endl;
    std::cout << "argv[2] = other one" << std::endl;
    std::cout << std::endl;
    std::exit(1);
  }

  std::cout << "Start" << std::endl;

  std::cout << "TFile::Open() file=" << argv[1] << std::endl;
  StartWatch();
  auto tf1 = TFile::Open(argv[1],"READ");
  StopWatch();

  std::cout << "TFile::Open() file=" << argv[2] << std::endl;
  StartWatch();
  auto tf2 = TFile::Open(argv[2],"READ");
  StopWatch();


 // KEY: TMatrixT<double> full_covariance_EXAMPLE1;1;
 // KEY: TMatrixT<double> frac_covariance_EXAMPLE1;1;
 // KEY: TMatrixT<double> full_correlation_EXAMPLE1;1;

  std::cout << std::endl;
  std::cout << "Read 3 matrix from file=" << tf1->GetName() << std::endl;
  StartWatch();
  auto full_covariance1  = (TMatrixD*)tf1->Get("full_covariance_EXAMPLE1");
  auto frac_covariance1  = (TMatrixD*)tf1->Get("frac_covariance_EXAMPLE1");
  auto full_correlation1 = (TMatrixD*)tf1->Get("full_correlation_EXAMPLE1");
  StopWatch();

  std::cout << std::endl;
  std::cout << "Read 3 matrix from file=" << tf2->GetName() << std::endl;
  StartWatch();
  auto full_covariance2  = (TMatrixD*)tf2->Get("full_covariance_EXAMPLE1");
  auto frac_covariance2  = (TMatrixD*)tf2->Get("frac_covariance_EXAMPLE1");
  auto full_correlation2 = (TMatrixD*)tf2->Get("full_correlation_EXAMPLE1");
  StopWatch();

  std::cout << std::endl;
  std::cout << "Compare covariance" << std::endl;
  Compare(full_covariance1,full_covariance2);

  std::cout << std::endl;
  std::cout << "Compare fractional" << std::endl;
  Compare(frac_covariance1,frac_covariance2);

  std::cout << std::endl;
  std::cout << "Compare correlation" << std::endl;
  Compare(full_correlation1,full_correlation2);

  global_watch.Stop();
  std::cout << "End t=(" << global_watch.CpuTime() << "," << global_watch.RealTime() << ")" << std::endl;
  return 0;
}
