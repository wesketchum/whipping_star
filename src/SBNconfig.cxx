#include "SBNconfig.h"
using namespace sbn;

SBNconfig::SBNconfig(std::vector<std::string> modein, std::vector<std::string> detin, std::vector<std::string> chanin, std::vector<std::vector<std::string>> subchanin, std::vector<std::vector<double>> binin){

	isVerbose = true;

	num_detectors = detin.size();
	num_channels = chanin.size();
	num_modes = modein.size();

	if(subchanin.size() != chanin.size()){
		std::cout<<"SUBCHAN.size() != chanin.size()"<<std::endl;
		exit(EXIT_FAILURE);
	}

	for(auto sb: subchanin){
		num_subchannels.push_back( sb.size());
	}

	for(auto bn: binin){
		num_bins.push_back(bn.size()-1);
	}

	this->calcTotalBins();

	xmlname = "NULL";	


	//the xml names are the way we track which channels and subchannels we want to use later
	mode_names = modein; 			
	detector_names = detin; 		
	channel_names = chanin; 		
	subchannel_names = subchanin; 

	for(auto c: chanin){
		channel_bool.push_back(true);
	}
	for(auto m: modein){
		mode_bool.push_back(true);
	}
	for(auto d: detin){
		detector_bool.push_back(true);
	}
	for(auto c: chanin){
		std::vector<bool> tml;
		for(int i=0; i< num_subchannels.size(); i++){
			tml.push_back(true);
		}
		subchannel_bool.push_back(tml);
	}

	//self explanatory
	bin_edges = binin;

	for(auto binedge: bin_edges){
		std::vector<double> binwidth;
		for(int b = 0; b<binedge.size()-1; b++){
			binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
		}
		bin_widths.push_back(binwidth);	
	}

	// The order is IMPORTANT its the same as defined in xml
	for(auto m: mode_names){
		for(auto d: detector_names){
			for(int c = 0; c< num_channels; c++){
				for(auto sb: subchannel_names.at(c)){
					std::string tmp = m +"_"+ d+"_"+channel_names.at(c)+"_"+sb;
					fullnames.push_back(tmp);
				}
			}
		}
	}


	for(int i=0; i< num_bins_total; i++){
		used_bins.push_back(i);
	}


};



SBNconfig::SBNconfig(std::string whichxml): xmlname(whichxml) {
	//standard constructor given an xml.
	//Using a very simple xml format that I directly coppied from an old project.

	isVerbose = true;
	has_oscillation_patterns = false;


	//max subchannels 100?
	subchannel_names.resize(100);
	subchannel_bool.resize(100);
	subchannel_osc_patterns.resize(100);
	char *end;

	//Setup TiXml documents
	TiXmlDocument doc( whichxml.c_str() );
	bool loadOkay = doc.LoadFile();
	if(loadOkay){
		std::cout<<"SBNconfig::SBNconfig || Loaded "<<whichxml<<std::endl;
	}else{
		std::cerr<<"ERROR: SBNonfig::SBNconfig || Failed to load "<<whichxml<<std::endl;
		exit(EXIT_FAILURE);
	}
	TiXmlHandle hDoc(&doc);


	// we have Modes, Detectors, Channels, Covariance matricies, MC multisim data, oscillation pattern matching
	TiXmlElement *pMode, *pDet, *pChan, *pCov, *pMC, *pData;


	//Grab the first element. Note very little error checking here! make sure they exist.
	pMode = doc.FirstChildElement("mode");
	pDet =  doc.FirstChildElement("detector");
	pChan = doc.FirstChildElement("channel");
	pCov  = doc.FirstChildElement("covariance");
	pMC   = doc.FirstChildElement("MCevents");
	pData   = doc.FirstChildElement("data");





	//Where is the "data" folder that keeps pre-computed spectra and rootfiles
	//Will eventuall just configure this in CMake
	while(pData){
		data_path = pData->Attribute("path");
		pData = pData->NextSiblingElement("data");
		std::cout<<"SBNconfig::SBnconfig || data path loaded as: "<<data_path<<std::endl;
	}


	// Where is the covariance matrix you want to use, and whats its name in the root file.
	while(pCov){
		correlation_matrix_rootfile = data_path + pCov->Attribute("file");
		correlation_matrix_name = pCov->Attribute("name");
		pCov = pCov->NextSiblingElement("covariance");
	}


	// What modes are we running in (e.g nu, nu bar, horn current=XXvolts....) Can have as many as we want
	while(pMode)
	{
		//std::cout<<"Mode: "<<pMode->Attribute("name")<<" "<<pMode->Attribute("use")<<std::endl;
		mode_names.push_back(pMode->Attribute("name"));	
		mode_bool.push_back(strtod(pMode->Attribute("use"),&end));	

		pMode = pMode->NextSiblingElement("mode");
		std::cout<<"SBNconfig::SBnconfig || loading mode: "<<mode_names.back()<<" with use_bool "<<mode_bool.back()<<std::endl;

	}


	// How many detectors do we want! 
	pDet = doc.FirstChildElement("detector");
	while(pDet)
	{
		//std::cout<<"Detector: "<<pDet->Attribute("name")<<" "<<pDet->Attribute("use")<<std::endl;
		detector_names.push_back(pDet->Attribute("name"));
		detector_bool.push_back(strtod(pDet->Attribute("use"),&end));
		pDet = pDet->NextSiblingElement("detector");	
		std::cout<<"SBNconfig::SBnconfig || loading detector: "<<detector_names.back()<<" with use_bool "<<detector_bool.back()<<std::endl;
	}

	//How many channels do we want! At the moment each detector must have all channels
	int nchan = 0;
	while(pChan)
	{

		// Read in how many bins this channel uses
		channel_names.push_back(pChan->Attribute("name"));
		channel_bool.push_back(strtod(pChan->Attribute("use"),&end));
		num_bins.push_back(strtod(pChan->Attribute("numbins"), &end));		
		std::cout<<"SBNconfig::SBNconfig || Loading Channel : "<<channel_names.back()<<" with use_bool: "<<channel_bool.back()<<std::endl;
		// What are the bin edges and bin widths (bin widths just calculated from edges now)
		TiXmlElement *pBin = pChan->FirstChildElement("bins");
		std::stringstream iss(pBin->Attribute("edges"));
		//std::stringstream pss(pBin->Attribute("widths"));

		double number;
		std::vector<double> binedge;
		std::vector<double> binwidth;
		while ( iss >> number ) binedge.push_back( number );
		//while ( pss >> number ) binwidth.push_back( number );
		for(int b = 0; b<binedge.size()-1; b++){
			binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
		}			
		bin_edges.push_back(binedge);
		bin_widths.push_back(binwidth);

		// Now loop over all this channels subchanels. Not the names must be UNIQUE!!
		TiXmlElement *pSubChan;

		pSubChan = pChan->FirstChildElement("subchannel");
		int nsubchan=0;
		while(pSubChan){
			//std::cout<<"Subchannel: "<<pSubChan->Attribute("name")<<" use: "<<pSubChan->Attribute("use")<<" osc: "<<pSubChan->Attribute("osc")<<std::endl;
			subchannel_names[nchan].push_back(pSubChan->Attribute("name"));
			subchannel_bool[nchan].push_back(strtod(pSubChan->Attribute("use"),&end));
			//0 means dont oscillate, 11 means electron disapearance, -11 means antielectron dis..etc..
			if(pSubChan->Attribute("osc"))
			{
				has_oscillation_patterns = true;
			}

			subchannel_osc_patterns.at(nchan).push_back(strtod(pSubChan->Attribute("osc"), &end));

			std::cout<<"--> Subchannel: "<<subchannel_names.at(nchan).back()<<" with use_bool "<<subchannel_bool.at(nchan).back()<<" and osc_pattern "<<subchannel_osc_patterns.at(nchan).back()<<std::endl;

			nsubchan++;
			pSubChan = pSubChan->NextSiblingElement("subchannel");	
		}
		num_subchannels.push_back(nsubchan);




		nchan++;
		pChan = pChan->NextSiblingElement("channel");	
	}

	// if wea re creating a covariance matrix using a ntuple and weights, here is the info
	if(pMC){
		std::cout<<"SBNcongig::SBNconfig || Loading a MC config. This is quite depreciated here."<<std::endl;
		while(pMC)
		{	
			pot.push_back(strtof(pMC->Attribute("pot"),&end));
			pot_scaling.push_back(strtof(pMC->Attribute("potscale"),&end));
			num_multisim.push_back(strtod(pMC->Attribute("multisim"),&end));
			multisim_name.push_back(pMC->Attribute("name"));
			multisim_file.push_back(pMC->Attribute("filename"));

			TiXmlElement *pParams = pMC->FirstChildElement("parameters");

			std::stringstream sss(pParams->Attribute("names"));

			std::vector<std::string> vstring;
			std::string nam;
			while ( sss >> nam) vstring.push_back( nam );

			parameter_names.push_back(vstring);



			TiXmlElement *pBranchT;
			pBranchT = pMC->FirstChildElement("btype");
			//			std::cout<<"Starting run over branch types"<<std::endl;
			std::vector<std::string> TEMP_branch_names_int;
			std::vector<std::string> TEMP_branch_names_int_array;
			std::vector<std::string> TEMP_branch_names_double;
			std::vector<std::string> TEMP_branch_names_double_array;
			std::vector<int> TEMP_branch_names_double_array_dimension;

			while(pBranchT){
				TiXmlElement *pBranch;
				pBranch = pBranchT->FirstChildElement("branch");

				while(pBranch){
					//std::cout<<pBranch->Attribute("name")<<" Type: "<<pBranchT->Attribute("type") <<std::endl;
					if( strtod(pBranchT->Attribute("type"),&end)  == 0){
						//std::cout<<pBranch->Attribute("name")<<" num_multisim"<<num_multisim<<" int"<<std::endl;
						TEMP_branch_names_int.push_back(pBranch->Attribute("name"));
					}
					else if (strtod(pBranchT->Attribute("type"),&end)  == 1){
						//std::cout<<pBranch->Attribute("name")<<" num_multisim"<<num_multisim<<" double"<<std::endl;
						TEMP_branch_names_double.push_back(pBranch->Attribute("name"));
					}
					else if (strtod(pBranchT->Attribute("type"),&end)  == 2){
						TEMP_branch_names_int_array.push_back(pBranch->Attribute("name"));
					}
					else if (strtod(pBranchT->Attribute("type"),&end)  == 3){
						//std::cout<<pBranch->Attribute("name")<<" num_multisim"<<num_multisim<<" double"<<std::endl;
						TEMP_branch_names_double_array.push_back(pBranch->Attribute("name"));
						//std::cout<<"Hi: "<<strtod(pBranch->Attribute("dimension"),&end)<<std::endl; 
						TEMP_branch_names_double_array_dimension.push_back(strtod(pBranch->Attribute("dimension"),&end));
					}
					pBranch = pBranch->NextSiblingElement("branch");	
				}
				pBranchT= pBranchT->NextSiblingElement("btype");
			}
			branch_names_double.push_back(TEMP_branch_names_double);
			branch_names_double_array.push_back(TEMP_branch_names_double_array);
			branch_names_double_array_dimension.push_back(TEMP_branch_names_double_array_dimension);

			branch_names_int_array.push_back(TEMP_branch_names_int_array);
			branch_names_int.push_back(TEMP_branch_names_int);
			pMC=pMC->NextSiblingElement("MCevents");
		}
	}





	// so num_channels here is number of TOTAL channels in xml.
	num_channels = channel_names.size();
	num_modes = mode_names.size();
	num_detectors  = detector_names.size();

	//Calculate bin_widths from bin_edges



	// here we run through every combination, and make note when (and where binwise) all the subchannels that are turned on are.
	std::string tempn;
	int indexcount = 0;
	for(int im = 0; im < num_modes; im++){

		for(int id =0; id < num_detectors; id++){

			for(int ic = 0; ic < num_channels; ic++){

				for(int sc = 0; sc < num_subchannels[ic]; sc++){

					tempn = mode_names[im] +"_" +detector_names[id]+"_"+channel_names[ic]+"_"+subchannel_names[ic][sc];

					// This is where you choose NOT to use some fields	
					if(mode_bool[im] && detector_bool[id] && channel_bool[ic] && subchannel_bool[ic][sc]){					
						fullnames.push_back(tempn);

						for(int k = indexcount; k < indexcount+num_bins[ic]; k++){
							used_bins.push_back(k);

						}
					}

					std::vector<int> tvec = {indexcount, indexcount+num_bins[ic]-1}; 

					mapIndex[tempn] = 	tvec;				
					indexcount = indexcount + num_bins[ic];


				}	
			}
		}
	}

	//Actually want num_XXX to contain only activated things
	//For here on down everything is derivable, above is just until I actually get config working.
	num_modes = 0;
	for(bool i:mode_bool){	if(i) num_modes++;	}

	num_detectors = 0;
	for(bool i:detector_bool){	if(i) num_detectors++;	}

	num_channels = 0;
	for(bool i:channel_bool){	if(i) num_channels++;	}

	for(int i=0; i< num_channels; i++){
		num_subchannels[i] = 0;
		for(bool j: subchannel_bool[i]){ if(j) num_subchannels[i]++;}
	}


	this->calcTotalBins();

}//end constructor


int SBNconfig::calcTotalBins(){

	// These variables are important
	// They show how big each mode block and decector block are, for any given number of channels/subchannels
	// both before and after compression!
	num_bins_detector_block = 0;
	num_bins_detector_block_compressed = 0;
	for(int i =0; i< num_channels; i++){
		num_bins_detector_block += num_subchannels[i]*num_bins[i];
		num_bins_detector_block_compressed += num_bins[i];	
	}		

	num_bins_mode_block = num_bins_detector_block*num_detectors;
	num_bins_mode_block_compressed = num_bins_detector_block_compressed*num_detectors;

	num_bins_total = num_modes*num_bins_mode_block;
	num_bins_total_compressed = num_modes*num_bins_mode_block_compressed;

	return 0;
}


