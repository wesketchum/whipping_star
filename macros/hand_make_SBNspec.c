void hand_make_SBNspec(){

TFile * f = new TFile ("unit1c.root","recreate");
f->cd();

std::vector<double> bins = {100,200};
TH1D *nu_SBND_nue_intrinsic = new TH1D("nu_SBND_nue_intrinsic","nu_SBND_nue_intrinsic",1,&bins[0]);
TH1D *nu_SBND_nue_leesignal = new TH1D("nu_SBND_nue_leesignal","nu_SBND_nue_leesignal",1,&bins[0]);

TH1D *nu_uBooNE_nue_intrinsic = new TH1D("nu_uBooNE_nue_intrinsic","nu_uBooNE_nue_intrinsic",1,&bins[0]);
TH1D *nu_uBooNE_nue_leesignal = new TH1D("nu_uBooNE_nue_leesignal","nu_uBooNE_nue_leesignal",1,&bins[0]);

TH1D *nu_ICARUS_nue_intrinsic = new TH1D("nu_ICARUS_nue_intrinsic","nu_ICARUS_nue_intrinsic",1,&bins[0]);
TH1D *nu_ICARUS_nue_leesignal = new TH1D("nu_ICARUS_nue_leesignal","nu_ICARUS_nue_leesignal",1,&bins[0]);

nu_SBND_nue_intrinsic->SetBinContent(1,22.1);
nu_uBooNE_nue_intrinsic->SetBinContent(1,15.9);
nu_ICARUS_nue_intrinsic->SetBinContent(1,11.4);

nu_SBND_nue_leesignal->SetBinContent(1,5.2);
nu_uBooNE_nue_leesignal->SetBinContent(1,4.1);
nu_ICARUS_nue_leesignal->SetBinContent(1,0.8);

nu_SBND_nue_intrinsic->Write();
nu_uBooNE_nue_intrinsic->Write();
nu_ICARUS_nue_intrinsic->Write();
nu_SBND_nue_leesignal->Write();
nu_uBooNE_nue_leesignal->Write();
nu_ICARUS_nue_leesignal->Write();


f->Close();



/*
TFile * f = new TFile ("unit1b.root","recreate");
f->cd();

std::vector<double> bins = {100,200};
TH1D *nu_uBooNE_nue_intrinsic = new TH1D("nu_uBooNE_nue_intrinsic","nu_uBooNE_nue_intrinsic",1,&bins[0]);
TH1D *nu_uBooNE_nue_leesignal = new TH1D("nu_uBooNE_nue_leesignal","nu_uBooNE_nue_leesignal",1,&bins[0]);

TH1D *nu_uBooNE_numu_intrinsic = new TH1D("nu_uBooNE_numu_intrinsic","nu_uBooNE_numu_intrinsic",1,&bins[0]);
TH1D *nu_uBooNE_numu_leesignal = new TH1D("nu_uBooNE_numu_leesignal","nu_uBooNE_numu_leesignal",1,&bins[0]);

TH1D *nu_uBooNE_nutau_intrinsic = new TH1D("nu_uBooNE_nutau_intrinsic","nu_uBooNE_nutau_intrinsic",1,&bins[0]);
TH1D *nu_uBooNE_nutau_leesignal = new TH1D("nu_uBooNE_nutau_leesignal","nu_uBooNE_nutau_leesignal",1,&bins[0]);


nu_uBooNE_nue_intrinsic->SetBinContent(1,22.1);
nu_uBooNE_numu_intrinsic->SetBinContent(1,15.9);
nu_uBooNE_nutau_intrinsic->SetBinContent(1,11.4);

nu_uBooNE_nue_leesignal->SetBinContent(1,5.2);
nu_uBooNE_numu_leesignal->SetBinContent(1,4.1);
nu_uBooNE_nutau_leesignal->SetBinContent(1,0.8);

nu_uBooNE_nue_intrinsic->Write();
nu_uBooNE_nue_leesignal->Write();
nu_uBooNE_numu_intrinsic->Write();
nu_uBooNE_numu_leesignal->Write();
nu_uBooNE_nutau_intrinsic->Write();
nu_uBooNE_nutau_leesignal->Write();



f->Close();
*/


/*
TFile * f = new TFile ("unit1a.root","recreate");
f->cd();

std::vector<double> bins = {100,200,300,400};
TH1D *nu_uBooNE_nue_intrinsic = new TH1D("nu_uBooNE_nue_intrinsic","nu_uBooNE_nue_intrinsic",3,&bins[0]);
TH1D *nu_uBooNE_nue_leesignal = new TH1D("nu_uBooNE_nue_leesignal","nu_uBooNE_nue_leesignal",3,&bins[0]);

nu_uBooNE_nue_intrinsic->SetBinContent(1,22.1);
nu_uBooNE_nue_intrinsic->SetBinContent(2,15.9);
nu_uBooNE_nue_intrinsic->SetBinContent(3,11.4);

nu_uBooNE_nue_leesignal->SetBinContent(1,5.2);
nu_uBooNE_nue_leesignal->SetBinContent(2,4.1);
nu_uBooNE_nue_leesignal->SetBinContent(3,0.8);

nu_uBooNE_nue_intrinsic->Write();
nu_uBooNE_nue_leesignal->Write();

f->Close();
*/
}
