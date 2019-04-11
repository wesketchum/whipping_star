void make_lee_friend_example_macro(){

    //Get the signal weights which are in the form of a TH1D
    TFile *flee = new TFile("/uboone/data/users/markross/SBNfit_example_data/LEE_combined_weights_wes_combined.root","read");
    TH1D * thlee = (TH1D*)flee->Get("LEE_combined_nue_CC_weights");

    //Get the nue's that we want to rescale
    TFile *fnue = new TFile("../1e1p.root","read");
    TTree *tnue = (TTree*)fnue->Get("events");
    //We would usually scale by true neutrino energy, however these simple test files do not have true neutrino energy so lets just use ereco for this example
    double energy;
    tnue->SetBranchAddress("ereco",&energy);

    //And create a file to store the outputted friend tree!
    TFile *ffriend = new TFile("/uboone/data/users/markross/SBNfit_example_data/lee_signal_friend_tree.root","recreate");
    TTree *tfriend = new TTree("lee_signal_weights","lee_signal_weights");

    //Create a branch in that tree to store the weights
    double lee_weights;
    tfriend->Branch("lee_weights",&lee_weights,"lee_weights/D");


    //Now we will loop over all events in the input nue file, getting the appropiate weight, and saving to the friend tree
    std::cout<<"Beginning to loop over "<<tnue->GetEntries()<<" entries in nue file"<<std::endl;
    for(int i=0; i<tnue->GetEntries(); i++){
        tnue->GetEntry(i);

        //Find which bin of LEE this entry lives in
        //NOTE everything here is in MeV
        int which_bin = thlee->GetXaxis()->FindBin(energy);
     
        if(i%100==0)std::cout<<"Nue energy: "<<energy<<" is bin#: "<<which_bin<<" ["<<thlee->GetBinLowEdge(which_bin)<<","<<thlee->GetBinLowEdge(which_bin)+thlee->GetBinWidth(which_bin)<<"] with weight: "<<thlee->GetBinContent(which_bin)<<std::endl;

        lee_weights = (double)thlee->GetBinContent(which_bin);

        //Fill once per entry in the nue ttree so that it matches up perfectly
        tfriend->Fill();
    }
    

    //And write out the 
    ffriend->cd();
    tfriend->Write();

    ffriend->Close();
    fnue->Close();
    flee->Close();
}


