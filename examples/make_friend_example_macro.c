void make_friend_example_macro(){



//Lets say you have a file, 1e1p.root containing a TTree called events. 
TFile *fin = new TFile("1e1p.root","read");
TTree *tin = (TTree*)fin->Get("events");

//We want to reweight this by some function of a branch, say ereco
double ereco;
tin->SetBranchAddress("ereco",&ereco);


//We want to store the new weights in another file and TTtee that we can easily assign as a FriendTree to the original
TFile *fout = new TFile("1e1p_weights.root","recreate");
fout->cd();
TTree *tout = new TTree("signalweights","signalweights");

//Lets create the weight that we will be filling
double wei;
tout->Branch("leeweight", &wei,"leeweight/D");


//Now we loop over all events in the input file, and calculate a weight based on our branch 
for(int i=0; i< tin->GetEntries(); ++i){
    tin->GetEntry(i);

    //This is just a nonsence quick test, you would put any weights appropiate here
    if(ereco > 600){ wei = 0.0;
    }else if(ereco > 500){
        wei = 1.0;
    }else if(ereco > 400){
        wei = 2.0;
    }else  if(ereco > 300){
        wei = 4.0;
    }

    //and fill the new tree. Because its being filled once for every event in the input TTree in unison, we know the entries match up and we know it will work well as a TTree.
    tout->Fill();
}

tout->Write();
fout->Close();


//The weights can then be easily accessed by assigning signalweights as a FriendTree via
// tin->AddFriend("signalweights","1e1p_weights.root");
// tin->Draw("ereco", "signalweights.leeweight");


return;
}
