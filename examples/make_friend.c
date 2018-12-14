void make_friend(){

TFile *fin = new TFile("1e1p.root","read");
TTree *tin = (TTree*)fin->Get("events");

double ereco;
tin->SetBranchAddress("ereco",&ereco);

TFile *fout = new TFile("1e1p_weights.root","recreate");
fout->cd();
TTree *tout = new TTree("signalweights","signalweights");

double wei;
tout->Branch("leeweight",&wei,"leeweight/D");

for(int i=0; i< tin->GetEntries(); ++i){
    tin->GetEntry(i);


    if(ereco > 600){ wei = 0.0;
    }else if(ereco > 500){
        wei = 1.0;
    }else if(ereco > 400){
        wei = 2.0;
    }else  if(ereco > 300){
        wei = 4.0;
    }

    tout->Fill();
}

tout->Write();
fout->Close();



return;
}
