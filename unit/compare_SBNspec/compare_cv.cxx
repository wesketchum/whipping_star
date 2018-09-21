#include <iostream>
#include <vector>
#include <sstream>


#include "TFile.h"
#include "TH1D.h"
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

void Compare(const TH1D* h1, const TH1D* h2) {

  int nbins = h1->GetNbinsX();

  StartWatch();
  for(int i=1; i<=nbins; ++i) {
    if ((*h1).GetBinContent(i) != (*h2).GetBinContent(i)) { 
      std::stringstream ss; 
      ss << "Difference" << std::endl;
      ss << "@(i)=("<<i<<") h1="<<(*h1).GetBinContent(i)<<" h2="<<(*h2).GetBinContent(i)<<std::endl;
      std::cout << ss.str() << std::endl;
      std::exit(1);
    }
  } 
  StopWatch();
  
}

int main(int argc, char** argv) {
  global_watch.Reset();
  global_watch.Start();

  if (argc != 3) {
    std::cout << std::endl;
    std::cout << "argv[1] = EXAMPLE1.SBNspec.root" << std::endl;
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

 // KEY: TH1D     nu_uBooNE_nue_intrinsic;1;
 // KEY: TH1D     nu_uBooNE_nue_leesignal;1;
 // KEY: TH1D     nu_uBooNE_numu_intrinsic;1;
 // KEY: TH1D     nu_uBooNE_ccpi0_intrinsic;1;

  std::cout << std::endl;
  std::cout << "Read 4 hist from file=" << tf1->GetName() << std::endl;
  StartWatch();
  auto nue_intrinsic1   = (TH1D*)tf1->Get("nu_uBooNE_nue_intrinsic");
  auto nue_leesignal1   = (TH1D*)tf1->Get("nu_uBooNE_nue_leesignal");
  auto numu_intrinsic1  = (TH1D*)tf1->Get("nu_uBooNE_numu_intrinsic");
  auto ccpi0_intrinsic1 = (TH1D*)tf1->Get("nu_uBooNE_ccpi0_intrinsic");
  StopWatch();

  std::cout << std::endl;
  std::cout << "Read 4 hist from file=" << tf2->GetName() << std::endl;
  StartWatch();
  auto nue_intrinsic2   = (TH1D*)tf2->Get("nu_uBooNE_nue_intrinsic");
  auto nue_leesignal2   = (TH1D*)tf2->Get("nu_uBooNE_nue_leesignal");
  auto numu_intrinsic2  = (TH1D*)tf2->Get("nu_uBooNE_numu_intrinsic");
  auto ccpi0_intrinsic2 = (TH1D*)tf2->Get("nu_uBooNE_ccpi0_intrinsic");
  StopWatch();

  std::cout << std::endl;
  std::cout << "Compare nue_intrinsic" << std::endl;
  Compare(nue_intrinsic1,nue_intrinsic2);

  std::cout << std::endl;
  std::cout << "Compare nue_leesignal" << std::endl;
  Compare(nue_leesignal1,nue_leesignal2);

  std::cout << std::endl;
  std::cout << "Compare numu_intrinsic" << std::endl;
  Compare(numu_intrinsic1,numu_intrinsic2);

  std::cout << std::endl;
  std::cout << "Compare ccpi0_intrinsic" << std::endl;
  Compare(ccpi0_intrinsic1,ccpi0_intrinsic2);

  global_watch.Stop();
  std::cout << std::endl;
  std::cout << "End t=(" << global_watch.CpuTime() << "," << global_watch.RealTime() << ")" << std::endl;
  return 0;
}
