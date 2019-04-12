int macro_convert_to_useful(){
    //Quick macro to convert the .txt's from the MiniBooNE data release into a useful SBNfit acceptable format
    //Will do 3 things.
    // (1) Transcribe the ntuple of fullosc into a root file.
    // (2) Generate nu_mu and nu_e "fake" spectra as drawn from their final spectra, litle harder
    // (3) Remove the statistical part of the fractional covariance and save as a TMatrixT<double>. Need to do (2) first.



    // (1) -------------
    // miniboone_numunuefullosc_ntuple.txt  (reconstructed neutrino energy (MeV) , true neutrino energy (MeV), neutrino baseline (cm) , and event weight)

    TFile * f1out= new TFile("MiniBooNE_fullosc_2018.root","recreate");
    f1out->cd();
    TTree * t1out= new TTree("events","events");

    double m_ereco=0; //in gev
    double m_etrue=0; //in gev
    double m_ltrue=0; //in km
    double m_wei=0;

    t1out->Branch("reco_energy",  &m_ereco);
    t1out->Branch("true_energy",  &m_etrue);
    t1out->Branch("true_baseline",&m_ltrue);
    t1out->Branch("pot_weight",   &m_wei);
    
    std::ifstream in;
    in.open("miniboone_numunuefullosc_ntuple.txt");

    int n=0;
    while (1) {
        n++;
        double Reco=0;
        double True=0; 
        double L = 0;
        double weight = 0;
        in >> Reco;
        in >> True;
        in >> L;
        in >> weight;   

        if (!in.good()) break;
        //std::cout<<n<<" R: "<<Reco/1000.0<<" T: "<<True/1000.0<<" L: "<<L/100000.0<<" weight: "<<weight<<std::endl;
        m_ereco = Reco/1000.0;
        m_etrue = True/1000.0;
        m_ltrue = L/100000.0;
        m_wei   = weight;
        t1out->Fill();

    }
    in.close();

    f1out->cd();
    t1out->Write();
    f1out->Close();

    //End of 1...


    // (2) -------------
    // miniboone_numu.txt and   miniboone_nuebgr_lowe.txt, values coppie directly
    // bins of NUE: 0.2 0.3 0.375 0.475 0.55  0.675  0.8  0.95  1.1 1.3  1.5 3.0
    // bins of NUMU: 0 0.500 0.700 0.900 1.100 1.300 1.500 1.700 1.900

    std::vector<double> nue_bins = {0.2, 0.3, 0.375, 0.475 ,0.55,  0.675 , 0.8 , 0.95 , 1.1, 1.3 , 1.5, 3.0};
    std::vector<double> nue_val = {361.002334,  216.002142,  239.436776,  127.517957,  179.035344,  133.901816,  139.020389,  113.446978,  81.204519,  98.603919 , 137.953204};
    double nue_tot = 0;
    for(double &v: nue_val) nue_tot+=v;


    std::vector<double> numu_bins = {0,0.500,0.700,0.900,1.100,1.300,1.500,1.700,1.900};
    std::vector<double> numu_val = {38564.217639,59339.405335,53069.519495,37171.337542,23002.153188,12423.361945,6012.845025,2801.295291};
     double numu_tot = 0;
    for(double &v: numu_val) numu_tot+=v;

    std::cout<<"Total number of events in nue_bkg is "<<nue_tot<<" and numu is "<<numu_tot<<std::endl;

    TH1D *h_nue = new TH1D("nue","nue",nue_val.size(),&nue_bins[0]);
    TH1D *h_numu = new TH1D("numu","numu",numu_val.size(),&numu_bins[0]);

    for(int i=0; i< nue_val.size(); i++){
        h_nue->SetBinContent(i+1, nue_val[i]);
    }

    for(int i=0; i< numu_val.size(); i++){
        h_numu->SetBinContent(i+1, numu_val[i]);
    }

 
    TH1D *t = new TH1D("test","test",numu_val.size(), &numu_bins[0]);


    TFile * f2out= new TFile("MiniBooNE_numu_2018.root","recreate");
    f2out->cd();
    TTree * t2out= new TTree("events","events");

    double m2_ereco=0; //in gev
    double m2_etrue=0; //in gev
    double m2_ltrue=0; //in km
    double m2_wei=0;

    t2out->Branch("reco_energy",  &m2_ereco);
    t2out->Branch("true_energy",  &m2_etrue);
    t2out->Branch("true_baseline",&m2_ltrue);
    t2out->Branch("pot_weight",   &m2_wei);
    
    int num_gen = 25000;


    for(int i=0; i< num_gen; i++){

        double tmp = h_numu->GetRandom();

        m2_ereco =  tmp;
        m2_etrue = m2_ereco;
        m2_ltrue = 0.5;
        m2_wei = numu_tot/(double)num_gen;
    
        t->Fill(tmp,numu_tot/(double)num_gen);
        t2out->Fill();
    }

    //A TEST
    /*
    TCanvas *c = new TCanvas();
    c->cd();
    h_numu->SetLineColor(kBlue);
    h_numu->Draw();

    t->SetLineStyle(9);
    t->SetLineColor(kRed);
    std::cout<<"Test integral "<<t->Integral()<<std::endl;
    t->Draw("same");
    c->Show();
    */ 
    
    f2out->cd();
    t2out->Write();
    f2out->Close();


    TFile * f3out= new TFile("MiniBooNE_nuebkg_2018.root","recreate");
    f3out->cd();
    TTree * t3out= new TTree("events","events");

    double m3_ereco=0; //in gev
    double m3_etrue=0; //in gev
    double m3_ltrue=0; //in km
    double m3_wei=0;

    t3out->Branch("reco_energy",  &m3_ereco);
    t3out->Branch("true_energy",  &m3_etrue);
    t3out->Branch("true_baseline",&m3_ltrue);
    t3out->Branch("pot_weight",   &m3_wei);
    

    for(int i=0; i< num_gen; i++){

        double tmp = h_nue->GetRandom();

        m3_ereco =  tmp;
        m3_etrue = m3_ereco;
        m3_ltrue = 0.5;
        m3_wei = numu_tot/(double)num_gen;
    
        t->Fill(tmp,numu_tot/(double)num_gen);
        t3out->Fill();
    }

    f3out->cd();
    t3out->Write();
    f3out->Close();


    // (3) Remove the statistical part of the fractional covariance and save as a TMatrixT<double>. Need to do (2) first.
    // miniboone_full_fractcovmatrix_nu_lowe.txt //fullosc/nue/mumu

    std::ifstream in2;
    in2.open("miniboone_full_fractcovmatrix_nu_lowe.txt");

    std::vector<std::vector<double>>v_mat;

    int tot_size = nue_val.size()*2+numu_val.size();
     n=0;
    while (1) {
        n++;

        std::vector<double> thisrow(tot_size,-99); 
        for(int i=0; i< tot_size; i++){
            in2 >> thisrow[i];
        }

        if (!in2.good()) break;
        //for(double &v:thisrow) std::cout<<v<<" ";
        //std::cout<<std::endl;
        v_mat.push_back(thisrow);

    }
    in2.close();

    TFile *f4 = new TFile("MiniBooNE_full_fractional_covariance_2018.root","recreate");
    f4->cd();
    TMatrixT<double> ans(tot_size,tot_size);
    for(int i=0; i< tot_size; i++){
    for(int j=0; j< tot_size; j++){

        ans(i,j) = v_mat[i][j];
        if(i==j && i >= nue_val.size() && i < nue_val.size()*2) ans(i,j) = ans(i,j) - nue_val[i];
        if(i==j && i >= nue_val.size()*2) ans(i,j) = ans(i,j) - numu_val[i];

    }
    }
   
    ans.Write("full_covariance",TObject::kWriteDelete);
    f4->Close();



    return 0;
}
