#pragma GCC optimize("O3","unroll-loops","inline")
#include "SBNosc.h"
#include <functional> 

using namespace sbn;

SBNosc::SBNosc(std::vector<TH1D> const & bghist, const char* xmldata) : SBNspec(bghist, xmldata) {
	working_model.zero();
	mass_step_size = 0.04;
	which_mode = BOTH_ONLY;
}
SBNosc::SBNosc(std::string name, const char* xmldata) : SBNspec(name, xmldata) {
	working_model.zero();
	mass_step_size = 0.04;
	which_mode = BOTH_ONLY;
}

SBNosc::SBNosc(std::string name, std::string whichxml) : SBNspec(name, whichxml) {
	working_model.zero();
	mass_step_size = 0.04;
	which_mode = BOTH_ONLY;
}

SBNosc::SBNosc(SBNspec & specin): SBNspec(specin){
	working_model.zero();
	which_mode = BOTH_ONLY;
}


SBNosc::SBNosc(std::string name, std::string whichxml, NeutrinoModel in) : SBNosc(name, whichxml) {
	LoadModel(in);

}

int SBNosc::LoadModel(NeutrinoModel in){
	working_model = in;
	calcMassSplittings();
	return 0;
}

void SBNosc::setModel(NeutrinoModel const & in){
    working_model = in;
    calcMassSplittings();
}

/*************************************************
 * for a given working_model.
 * Calculate how many mass splittings, and which "type" it is,
 *
 * **********************************************/

int SBNosc::calcMassSplittings(){
        mass_splittings.clear();

	if (working_model.numsterile == 1) {
            mass_splittings.push_back(std::make_pair(0.0, 41));
            return 0;
	}

        // NOTE: the pair's first item is use nowehere ...

	double fix41   = round(log10( working_model.dm41Sq)/mass_step_size)*mass_step_size;
	double fix51   = round(log10(    (working_model.dm51Sq))/mass_step_size)*mass_step_size;
	double fix61   = round(log10(    (working_model.dm61Sq))/mass_step_size)*mass_step_size;

	double round54 = round(log10(fabs(working_model.dm54Sq))/mass_step_size)*mass_step_size;
	double round64 = round(log10(fabs(working_model.dm64Sq))/mass_step_size)*mass_step_size;
	double round65 = round(log10(fabs(working_model.dm65Sq))/mass_step_size)*mass_step_size;

       	if (working_model.numsterile ==2)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (round54,54));
	}
	else if (working_model.numsterile ==3)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (fix61,61));

		mass_splittings.push_back( std::make_pair (round54,54));
		mass_splittings.push_back( std::make_pair (round64,64));
		mass_splittings.push_back( std::make_pair (round65,65));
	}

	return 0;
}




int SBNosc::OscillateThis(std::string tag){
	this->CalcFullVector();
	this->CollapseVector();

	calcMassSplittings();

	for(auto ms: mass_splittings){

			//this is wrong
			std::string name_sinsq = tag +"_SINSQ_dm_"+working_model.mass_tag+".SBNspec.root";
			std::string name_sin = tag +"_SIN_dm_"+working_model.mass_tag+".SBNspec.root";


			SBNspec single_frequency(name_sin , xmlname , false);
			SBNspec single_frequency_square(name_sinsq , xmlname ,false);


			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu =working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu = 0.0;
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;
						prob_muebar_sq = 0;
						break;

				}



			//std::cout<<"mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;


			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp );
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
				}
			}



			/*single_frequency.Scale("elike_fulloscnue", prob_mue);
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq);
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/

			this->Add(&single_frequency);
			this->Add(&single_frequency_square);


	}//Done looping over

	this->CalcFullVector();
	this->CollapseVector();


	return 0;
};

std::vector<double> SBNosc::Oscillate(std::string tag, double scale){

	std::vector<double> tmp = this->Oscillate(tag);
	for(auto & v: tmp){
		v=v*scale;
	}

	return tmp;
}

/*
std::vector<double> SBNosc::OscillateWithAmp(double amp, double amp_sq){

		this->CalcFullVector();
		this->CollapseVector();

	std::vector<double > temp = collapsed_vector;


	for(auto ms: mass_splittings){
			char namei[200];

			sPrintf(namei, "%sprecomp/SBN_SIN_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency(namei , xmlname, false);

			sPrintf(namei, "%sprecomp/SBN_SINSQ_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency_square(namei , xmlname, false);

			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 0:
							break;
						default:
							osc_amp = amp;
							osc_amp_sq = amp_sq;
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );

				}
			}


			single_frequency.CalcFullVector();
			single_frequency.CollapseVector();


			single_frequency_square.CalcFullVector();
			single_frequency_square.CollapseVector();

			for(int i=0;i<temp.size(); i++){
				temp[i] += single_frequency.collapsed_vector[i];
				temp[i] += single_frequency_square.collapsed_vector[i];
			}

	}//Done looping over

	return temp;
};
*/

std::vector<double> SBNosc::Oscillate(std::string tag){
    return this->Oscillate(tag,true);
}


void scaleSubBy(std::vector<double> & vin, size_t first, size_t last, double scaleby) {
  std::transform(vin.begin() + first, vin.begin() + last, vin.begin() +first, 
      std::bind(std::multiplies<double>(), std::placeholders::_1, scaleby));
}

std::vector<double> SBNosc::Oscillate(std::vector<double> const & sf_sinsq, std::vector<double> const & sf_sin) {

    std::vector<double> out = {0.0, sf_sin.size()};
    std::vector<double> v1, v2;
    for (auto ms: mass_splittings) {
        int which_dm = ms.second;
        double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

        // TODO this should be in a function.
        switch (which_mode) {
            case APP_ONLY: //Strictly nu_e app only
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DIS_ONLY: //Strictly nu_mu dis only
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                break;
            case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                break;
        }

        v1 = sf_sinsq;
        v2 = sf_sin;
        double osc_amp(0), osc_amp_sq(0);
        int osc_pattern(0);
        // Iterate over channels
        size_t offset(0);
        for (int i=0; i<num_channels; i++) {
            size_t nbins_chan = num_bins.at(i);
            auto const & thisPattern = subchannel_osc_patterns[i];//.at(j);
            for (int j=0; j<num_subchannels[i]; j++){
                osc_pattern = thisPattern[j];//_osc_patterns.at(i).at(j);
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<num_detectors;++d) {
                  size_t first  = d*num_bins_detector_block + offset;
                  scaleSubBy(v1, first, first + nbins_chan, osc_amp_sq);
                  scaleSubBy(v2, first, first + nbins_chan, osc_amp   );
                }
                offset +=nbins_chan;
            }
	}
        std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(out), std::plus<double>());

    } // Done looping over mass splittings
    return out;
};

Eigen::VectorXd SBNosc::Oscillate(Eigen::VectorXd const & sf_sinsq, Eigen::VectorXd const & sf_sin) {
    Eigen::VectorXd retVec(num_bins_total);
    retVec.setZero(); // !!!

    for (auto ms: mass_splittings) {
        int which_dm = ms.second;
        double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

        // TODO this should be in a function.
        //switch (which_mode) {
            //case APP_ONLY: //Strictly nu_e app only
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                //break;
            //case DIS_ONLY: //Strictly nu_mu dis only
                //prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                //break;
            //case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                //prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                //prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                //prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                //prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                //prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                //prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                //break;
            //case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                //prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                //prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                //prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                //prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                //prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                //break;
            //case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                //prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                //break;
        //}

        double osc_amp(0), osc_amp_sq(0);
        int osc_pattern(0);
        // Iterate over channels
        size_t offset(0);
        for (int i=0; i<num_channels; i++) {
            size_t nbins_chan = num_bins[i];
            auto const & thisPattern = subchannel_osc_patterns[i];//.at(j);
            for (int j=0; j<num_subchannels[i]; j++){
                osc_pattern = thisPattern[j];//_osc_patterns.at(i).at(j);
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<num_detectors;++d) {
                  size_t first  = d*num_bins_detector_block + offset;
                  retVec.segment(first, nbins_chan) += osc_amp   *  sf_sin.segment(first, nbins_chan);
                  retVec.segment(first, nbins_chan) += osc_amp_sq*sf_sinsq.segment(first, nbins_chan);
                  //Eigen::Map<Eigen::VectorXd>(sf_sin.data()   + first, nbins_chan,1) *= osc_amp;
                  //Eigen::Map<Eigen::VectorXd>(sf_sinsq.data() + first, nbins_chan,1) *= osc_amp_sq;
                }
                offset +=nbins_chan;
            }
	}
        //retVec += sf_sin;
        //retVec += sf_sinsq;

    } // Done looping over mass splittings
    return retVec;
};

// This is a stripped down version of the code that really only does the
// oscillation. To get a final spectrum, the core spectrum needs to be added.
std::vector<double> SBNosc::Oscillate(
          std::unordered_map <std::string, Eigen::VectorXd > const & sinsqmap,
          std::unordered_map <std::string, Eigen::VectorXd > const & sinmap) 
{
    Eigen::VectorXd retVec(num_bins_total);
    retVec.setZero();


    for (auto ms: mass_splittings) {
        int which_dm = ms.second;
        double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

        // TODO this should be in a function.
        switch (which_mode) {
            case APP_ONLY: //Strictly nu_e app only
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DIS_ONLY: //Strictly nu_mu dis only
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                break;
            case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                break;
        }

        Eigen::VectorXd sf_sinsq = sinsqmap.at(working_model.mass_tag);
        Eigen::VectorXd sf_sin   = sinmap.at(  working_model.mass_tag);

        double osc_amp(0), osc_amp_sq(0);
        int osc_pattern(0);
        // Iterate over channels
        size_t offset(0);
        for (int i=0; i<num_channels; i++) {
            size_t nbins_chan = num_bins.at(i);
            for (int j=0; j<num_subchannels.at(i); j++){
                osc_pattern = subchannel_osc_patterns.at(i).at(j);
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<num_detectors;++d) {
                  size_t first  = d*num_bins_detector_block + offset;
                  Eigen::Map<Eigen::VectorXd>(sf_sin.data()   + first, nbins_chan,1) *= osc_amp;
                  Eigen::Map<Eigen::VectorXd>(sf_sinsq.data() + first, nbins_chan,1) *= osc_amp_sq;
                }
                offset +=nbins_chan;
            }
	}
        retVec += sf_sin;
        retVec += sf_sinsq;

    } // Done looping over mass splittings
    std::vector<double> temp(retVec.data(), retVec.data() + retVec.rows() * retVec.cols());
    //for (int i=0; i<num_bins_total;++i) std::cerr << "EIG :" << i << " " << temp[i] << "\n";
    return temp;
};

// This is a stripped down version of the code that really only does the
// oscillation. To get a final spectrum, the core spectrum needs to be added.
std::vector<double> SBNosc::Oscillate(
          std::unordered_map <std::string, std::vector<double> > const & sinsqmap,
          std::unordered_map <std::string, std::vector<double> > const & sinmap) 
{
    Eigen::VectorXd retVec(num_bins_total);
    retVec.setZero();


    for (auto ms: mass_splittings) {
        int which_dm = ms.second;
        double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

        // TODO this should be in a function.
        switch (which_mode) {
            case APP_ONLY: //Strictly nu_e app only
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DIS_ONLY: //Strictly nu_mu dis only
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                break;
            case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                break;
        }

        std::vector<double> sf_sinsq = sinsqmap.at(working_model.mass_tag);
        std::vector<double> sf_sin   = sinmap.at(  working_model.mass_tag);

        double osc_amp(0), osc_amp_sq(0);
        int osc_pattern(0);
        // Iterate over channels
        size_t offset(0);
        for (int i=0; i<num_channels; i++) {
            size_t nbins_chan = num_bins.at(i);
            for (int j=0; j<num_subchannels.at(i); j++){
                osc_pattern = subchannel_osc_patterns.at(i).at(j);
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<num_detectors;++d) {
                  size_t first  = d*num_bins_detector_block + offset;
                  size_t last   = d*num_bins_detector_block + offset + nbins_chan;
                  scaleSubBy(sf_sinsq, first, last, osc_amp_sq);
                  scaleSubBy(sf_sin,   first, last, osc_amp   );
                }
                offset +=nbins_chan;
            }
	}
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sf(sf_sin.data(), num_bins_total, 1);
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sfsq(sf_sinsq.data(), num_bins_total, 1);
        retVec += _sf;
        retVec += _sfsq;

    } // Done looping over mass splittings
    std::vector<double> temp(retVec.data(), retVec.data() + retVec.rows() * retVec.cols());
    return temp;
};

// Use vector< vector<double>> instead
std::vector<double> SBNosc::Oscillate(std::string tag, bool return_compressed, const char * xmldata,
          std::unordered_map <std::string, std::vector<TH1D> > const & sinsqmap,
          std::unordered_map <std::string, std::vector<TH1D> > const & sinmap) 
{
    this->CalcFullVector();
    this->CollapseVector();
    calcMassSplittings();
    std::vector<double> temp;
    
    if (return_compressed)  temp = collapsed_vector;
    else temp = full_vector;

    int n = temp.size();
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > PRED(temp.data(), n, 1);
    //PRED.setZero(); // uncomment to prevent adding of the core spectrum

    for (auto ms: mass_splittings) {

        SBNspec single_frequency_square(sinsqmap.at(working_model.mass_tag) , xmldata , false);
        SBNspec single_frequency(       sinmap.at(  working_model.mass_tag) , xmldata , false);

        single_frequency.CalcFullVector();
        single_frequency_square.CalcFullVector();

        // This is madness, we only want two numbers, osc_amp and osc_amp_sq as function of dm.
        // they are used to scale (i.e. multiply) the spectra with

        int which_dm = ms.second;
        double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

        // TODO this should be in a function.
        switch (which_mode) {
            case APP_ONLY: //Strictly nu_e app only
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DIS_ONLY: //Strictly nu_mu dis only
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                break;
            case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                prob_mumu      = working_model.oscAmp( 2,  2, which_dm, 2);
                prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
                prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
                prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
                prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
                break;
            case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                prob_ee        = working_model.oscAmp( 1,  1, which_dm, 2);
                break;
        }

       double osc_amp(0), osc_amp_sq(0);
       int osc_pattern(0);
       for (int i=0; i<num_channels; i++){
           for (int j=0; j<num_subchannels.at(i); j++){
               osc_pattern = subchannel_osc_patterns.at(i).at(j);
               switch (osc_pattern){
                   case 11:
                       osc_amp_sq = prob_ee;
                       break;
                   case -11:
                       osc_amp_sq = prob_ee;
                       break;
                   case 22:
                       osc_amp_sq = prob_mumu;
                       break;
                   case -22:
                       osc_amp_sq = prob_mumu;
                       break;
                   case 21:
                       osc_amp    = prob_mue;
                       osc_amp_sq = prob_mue_sq;
                       break;
                   case -21:
                       osc_amp    = prob_muebar;
                       osc_amp_sq = prob_muebar_sq;
                       break;
                   case 0:
                   default:
                       break;
               }

               single_frequency.Scale(       channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp    ); // This calls TH1D on a subset of the histograms
               single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
           }
       }

       single_frequency.CalcFullVector();
       single_frequency.CollapseVector();

       single_frequency_square.CalcFullVector();
       single_frequency_square.CollapseVector();

       if (return_compressed) {

           Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sf(single_frequency.collapsed_vector.data(), n, 1);
           Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sfsq(single_frequency_square.collapsed_vector.data(), n, 1);
           PRED += _sf;
           PRED +=_sfsq;
       }
       else {
           Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sf(single_frequency.full_vector.data(), n, 1);
           Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > _sfsq(single_frequency_square.full_vector.data(), n, 1);
           PRED +=_sf;
           PRED +=_sfsq;
       }
    }//Done looping over
    std::vector<double> ret(PRED.data(), PRED.data() + PRED.rows() * PRED.cols());
    //for (int i=0; i<num_bins_total;++i) std::cerr << "OLD :" << i << " " << ret[i] << "\n";

    return ret;
};


std::vector<double> SBNosc::Oscillate(std::string tag, bool return_compressed, const char * xmldata) {
    this->CalcFullVector();
    this->CollapseVector();
    calcMassSplittings();
    std::vector<double> temp;
    
    if (return_compressed)  temp = collapsed_vector;
    else temp = full_vector;

    for (auto ms: mass_splittings) {

              std::string name_sinsq = tag +"_SINSQ_dm_"+working_model.mass_tag+".SBNspec.root";
              std::string name_sin = tag +"_SIN_dm_"+working_model.mass_tag+".SBNspec.root";

              // TODO replace with ctor that takes vector<TH1D>
              SBNspec single_frequency(name_sin , xmldata , false);
              SBNspec single_frequency_square(name_sinsq , xmldata ,false);

              single_frequency.CalcFullVector();
              single_frequency_square.CalcFullVector();

              double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

              int which_dm = ms.second;

              switch (which_mode)
                    {
                            case APP_ONLY: //Strictly nu_e app only
                                    prob_mumu =0;
                                    prob_ee   =0;
                                    prob_mue = working_model.oscAmp(2,1,which_dm,1);
                                    prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
                                    prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
                                    prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
                                    break;
                            case DIS_ONLY: //Strictly nu_mu dis only
                                    prob_mumu = working_model.oscAmp(2,2,which_dm,2);
                                    prob_ee = 0;
                                    prob_mue = 0;
                                    prob_mue_sq =0;
                                    prob_muebar =0;
                                    prob_muebar_sq =0;
                                    break;
                            case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
                                    prob_mumu = working_model.oscAmp(2,2,which_dm,2);
                                    prob_ee = working_model.oscAmp(1,1,which_dm,2);
                                    prob_mue = working_model.oscAmp(2,1,which_dm,1);
                                    prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
                                    prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
                                    prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
                                    break;
                            case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
                                    prob_mumu =working_model.oscAmp(2,2,which_dm,2);
                                    prob_ee = 0.0;
                                    prob_mue = working_model.oscAmp(2,1,which_dm,1);
                                    prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
                                    prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
                                    prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
                                    break;
                            case DISE_ONLY: // A strange version where nu_e can appear but not disapear
                                    prob_mumu = 0.0;
                                    prob_ee = working_model.oscAmp(1,1,which_dm,2);
                                    prob_mue = 0;
                                    prob_mue_sq = 0;
                                    prob_muebar = 0;
                                    prob_muebar_sq = 0;
                                    break;

                    }


			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
                                        int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
				}
			}

			single_frequency.CalcFullVector();
			single_frequency.CollapseVector();

			single_frequency_square.CalcFullVector();
			single_frequency_square.CollapseVector();



           if (return_compressed) { 
              for(int i=0;i<temp.size(); i++){
                temp[i] += single_frequency.collapsed_vector[i];
                temp[i] += single_frequency_square.collapsed_vector[i];
              }
           }
           else {
              for(int i=0;i<temp.size(); i++){
                temp[i] += single_frequency.full_vector[i];
                temp[i] += single_frequency_square.full_vector[i];
              }
           }
	}//Done looping over
	return temp;
};


std::vector<double> SBNosc::Oscillate(std::string tag, bool return_compressed){

		this->CalcFullVector();
		this->CollapseVector();

	calcMassSplittings();

	std::vector<double> temp;
    
    if(return_compressed){
        temp = collapsed_vector;
    }else {
        temp = full_vector;
    }

	for(auto ms: mass_splittings){

			std::string name_sinsq = tag +"_SINSQ_dm_"+working_model.mass_tag+".SBNspec.root";
			std::string name_sin = tag +"_SIN_dm_"+working_model.mass_tag+".SBNspec.root";

			SBNspec single_frequency(name_sin , xmlname , false);
			SBNspec single_frequency_square(name_sinsq , xmlname ,false);

			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

    		single_frequency.CalcFullVector();
			single_frequency_square.CalcFullVector();
            //single_frequency.PrintFullVector();
            //single_frequency_square.PrintFullVector();


			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu =working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu = 0.0;
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;
						prob_muebar_sq = 0;
						break;

				}


			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
				}
			}


	    	//std::cout<<"mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;
			/*
			single_frequency.Scale("elike_fulloscnue", prob_mue);
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq);
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/


			single_frequency.CalcFullVector();
			single_frequency.CollapseVector();

			single_frequency_square.CalcFullVector();
			single_frequency_square.CollapseVector();

           if(return_compressed){ 
			for(int i=0;i<temp.size(); i++){
				temp[i] += single_frequency.collapsed_vector[i];
				temp[i] += single_frequency_square.collapsed_vector[i];
			}

           }else{
        	for(int i=0;i<temp.size(); i++){
				temp[i] += single_frequency.full_vector[i];
				temp[i] += single_frequency_square.full_vector[i];
			}



           }

	}//Done looping over

	return temp;
};




int SBNosc::SetMode(int in){
	which_mode = in;

return in;
}

void SBNosc::SetAppMode(){
	SetMode(APP_ONLY);
}

void SBNosc::SetDisMode(){
	SetMode(DIS_ONLY);
}

void SBNosc::SetBothMode(){
	SetMode(BOTH_ONLY);
}

void SBNosc::SetWierdMode(){
	SetMode(WIERD_ONLY);
}

void SBNosc::SetDisEMode(){
	SetMode(DISE_ONLY);
}
