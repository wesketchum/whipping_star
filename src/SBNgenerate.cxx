#include "SBNgenerate.h"
//#include "MCEventWeight.h"

using namespace sbn;



SBNgenerate::SBNgenerate(std::string xmlname) {
    NeutrinoModel nullModel;
    SBNgenerate(xmlname, nullModel);
}

SBNgenerate::SBNgenerate(std::string xmlname, NeutrinoModel inModel ) : SBNconfig(xmlname), nu_model(inModel) {

    std::vector<double> FUDGE ={470.0, 470.0,600.0, 600.0, 600.0,100.0,100.0 };


    TRandom3 *rangen = new TRandom3(0);

    //	gROOT->ProcessLine("#include <map>");
    //	gROOT->ProcessLine("#include <vector>");
    //	gROOT->ProcessLine("#include <string>");

    //gSystem->Load("../src/libranch_weightsMapDict.so");

    //	std::string dict_location = "../../dict/AutoDict_map_string__vector_double____cxx.so";
    //	gSystem->Load(  (dict_location).c_str());

    //	gSystem->Load("/uboone/app/users/markrl/sbnfit/whipping_star/src/mdict_h.so");
    //
    //std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
    //

    std::map<std::string, int> parameter_sims;

    //Initialise the central value SBNspec.
    SBNspec tm(xmlname,-1,false);
    spec_central_value = tm;
    spec_osc_sin  = tm;
    spec_osc_sinsq  = tm;

    int num_files = montecarlo_file.size();

    montecarlo_additional_weight.resize(num_files,1.0);


    //Some POT counting

    for(auto &fn: montecarlo_file){
        files.push_back(new TFile(fn.c_str()));
    }


    for(int i=0; i<montecarlo_name.size(); i++){
        trees.push_back((TTree*)files.at(i)->Get(montecarlo_name.at(i).c_str()) );
    }

    for(int i=0; i<montecarlo_file.size(); i++){

        if( montecarlo_file_friend_treename_map.count(montecarlo_file.at(i))>0){
            for(int k=0; k< montecarlo_file_friend_treename_map.at(montecarlo_file.at(i)).size(); k++){

                std::string treefriendname = (montecarlo_file_friend_treename_map.at(montecarlo_file.at(i))).at(k);
                std::string treefriendfile = (montecarlo_file_friend_map.at(montecarlo_file.at(i))).at(k);

                std::cout<<"SBNmontecarlo::SBNmontecarlo\t|| Adding a friend tree  "<< treefriendfile<<" to file "<<montecarlo_file.at(i)<<std::endl;

                trees.at(i)->AddFriend( treefriendname.c_str()   ,  treefriendfile.c_str()   );
            }
        }
    }

    std::vector<int> nentries;
    for(auto &t: trees){
        nentries.push_back(t->GetEntries());
    }

    f_weights.resize(num_files,nullptr);



    for(int i=0; i< num_files; i++){

        double pot_scale = 1.0;
        if(montecarlo_pot[i]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[i];
        }

        montecarlo_scale[i] = montecarlo_scale[i]*pot_scale;




        //delete f_weights->at(i);	f_weights->at(i) = 0;
        trees[i]->SetBranchAddress("eventweights", &(f_weights[i]) );
        //delete f_weights->at(i);	f_weights->at(i) = 0;
        
        for(int k=0; k<branch_variables.at(i).size(); k++){
            std::cout<<"Setting Branch: "<<branch_variables.at(i).at(k)->name<<std::endl;
            trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->name.c_str(), branch_variables.at(i).at(k)->GetValue() );

            if(branch_variables.at(i).at(k)->GetOscillate()){
                std::cout<<"Setting true branch variables"<<std::endl;
                trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->true_param_name.c_str(), branch_variables.at(i).at(k)->GetTrueValue() );
                trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->true_L_name.c_str(), branch_variables.at(i).at(k)->GetTrueL() );
            }
        }

        if(montecarlo_additional_weight_bool[i]){
            //we have an additional weight we want to apply at run time, otherwise its just set at 1. 
            trees[i]->SetBranchAddress(montecarlo_additional_weight_names[i].c_str(), &montecarlo_additional_weight[i]); 
        }


    }
    std::string     bnbcorrection_str = "bnbcorrection_FluxHist";


    std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
    std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
    std::vector<double> base_vec (spec_central_value.num_bins_total,0.0);

    for(int j=0;j<num_files;j++){

        //delete f_weights->at(j);
        //f_weights->at(j)=0;

        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);
            std::map<std::string,std::vector<float>>* thisfWeight = f_weights[j];

            if(i%100==0) std::cout<<"SBNgenerate::SBNgenerate\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = montecarlo_additional_weight[j];
            global_weight = global_weight*montecarlo_scale[j];


            if(thisfWeight->count("bnbcorrection_FluxHist")>0){
                //global_weight = global_weight*thisfWeight->at("bnbcorrection_FluxHist").front();
            }


            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<"SBNgenerate::SBNgenerate\t|| ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }

            if( this->EventSelection(j) ){
                for(int t=0; t<branch_variables[j].size();t++){
                    //std::cout<<"Starting branch : "<<branch_variables.at(j).at(t)->name<<" "<<branch_variables.at(j).at(t)->associated_hist<<std::endl;
                    //Need the histogram index, the value, the global bin...

                    int ih = spec_central_value.map_hist.at(branch_variables[j][t]->associated_hist);
                    double reco_var = *(static_cast<double*>(branch_variables[j][t]->GetValue()));
                    int reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);


//                    reco_var = reco_var*1.031;

                    //std::cout<<ih<<" "<<reco_var<<" "<<reco_bin<<" JJ"<<std::endl;

                    //Find if this event should be oscillated
                    if(branch_variables[j][t]->GetOscillate()){
                        //Working
                        double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));
                        double true_L = *(static_cast<double*>(branch_variables[j][t]->GetTrueL()));

                        //first subtract off a random 50m
                        //true_L = true_L - rangen->Uniform(0,50.0);

                        if(ih==0 || ih ==1) true_L = true_L-10.0;

                        //WARNING need to change to km
                        true_L = true_L/1000.0;
                        //true_L = FUDGE[ih]/1000.0;
                        //
                        //std::cout<<ih<<" "<<spec_osc_sin.hist[ih].GetName()<<std::endl;

                        double osc_Probability_sin = nu_model.oscProbSin(true_var, true_L);
                        double osc_Probability_sinsq = nu_model.oscProbSinSq(true_var, true_L);

                        spec_osc_sinsq.hist[ih].Fill(reco_var, global_weight*osc_Probability_sinsq);
                        spec_osc_sin.hist[ih].Fill(reco_var, global_weight*osc_Probability_sin);
                        spec_central_value.hist[ih].Fill(reco_var,global_weight);
                        //std::cout<<"Reco: "<<reco_var<<" True: "<<true_var<<" L: "<<true_L<<" "<<osc_Probability_sin<<" "<<osc_Probability_sinsq<<" glob: "<<global_weight<<std::endl;
                    }else{
                        std::cout<<"Not Oscillated"<<std::endl;
                        spec_central_value.hist[ih].Fill(reco_var,global_weight);
                        spec_osc_sinsq.hist[ih].Fill(reco_var, global_weight);
                        spec_osc_sin.hist[ih].Fill(reco_var, global_weight);
                        //	std::cout<<reco_var<<" "<<std::endl;
                    }
                }
            }
        } //end of entry loop
    }//end of file loop


    /***************************************************************
     *		Now some clean-up and Writing
     * ************************************************************/

    for(auto &f: files){
        f->Close();
    }

}
/***************************************************************
 *		Some virtual functions for selection and histogram filling
 * ************************************************************/

int SBNgenerate::WritePrecomputedOscSpecs(std::string tag){

    std::cout<<"SBNGenerate::WritePrecomputedOscSpecs()\t\t||\t\tWriting out "<<tag<<"SINXX_dm_"<<nu_model.mass_tag<<std::endl;
    spec_osc_sinsq.WriteOut(tag+"_SINSQ_dm_"+nu_model.mass_tag);
    spec_osc_sin.WriteOut(tag+"_SIN_dm_"+nu_model.mass_tag);

    return 0;
}


bool SBNgenerate::EventSelection(int which_file){
    return true;
}

int SBNgenerate::FillHistograms(int file, int uni, double wei){
    return 0;
}
