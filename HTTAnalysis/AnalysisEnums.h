#ifndef RootAnalysis_AnalysisEnums_H
#define RootAnalysis_AnalysisEnums_H

#include <map>
#include <vector>
#include <iostream>

namespace HTTAnalysis {

///Copy from DataFormats/TauReco/interface/PFTauDecayMode.h
enum hadronicTauDecayModes
{
        tauDecay1ChargedPion0PiZero,
        tauDecay1ChargedPion1PiZero, // rho (770 MeV) mediated)
        tauDecay1ChargedPion2PiZero, // a1  (1.2 GeV) mediated
        tauDecay1ChargedPion3PiZero, // contaminated or unmerged photo
        tauDecay1ChargedPion4PiZero, // contaminated or unmerged photo
        tauDecay2ChargedPion0PiZero, // extra track or un-recod track
        tauDecay2ChargedPion1PiZero, // extra track or un-recod track
        tauDecay2ChargedPion2PiZero, // extra track or un-recod track
        tauDecay2ChargedPion3PiZero, // extra track or un-recod track
        tauDecay2ChargedPion4PiZero, // extra track or un-recod track
        tauDecay3ChargedPion0PiZero, // a1  (1.2 GeV) mediated
        tauDecay3ChargedPion1PiZero, // a1  (1.2 GeV) mediated
        tauDecay3ChargedPion2PiZero, // a1  (1.2 GeV) mediated
        tauDecay3ChargedPion3PiZero, // a1  (1.2 GeV) mediated
        tauDecay3ChargedPion4PiZero, // a1  (1.2 GeV) mediated
        tauDecaysElectron,
        tauDecayMuon,
        tauDecayOther             // catch-all
};

enum eventCategories {jet0_low, jet0_high,
                    jet1_low, jet1_high,
                    vbf_low, vbf_high,
                    jet0, boosted, vbf,
                    wjets_jet0, wjets_boosted, wjets_vbf,
                    wjets_qcd_jet0, wjets_qcd_boosted, wjets_qcd_vbf,
                    qcd_jet0, qcd_boosted, qcd_vbf,
                    qcd_ss_jet0, qcd_ss_boosted, qcd_ss_vbf,
                    ss_jet0, ss_boosted, ss_vbf,
                    antiIso_jet0, antiIso_boosted, antiIso_vbf,
                    mu_pi, mu_rho,  pi_pi, pi_rho, rho_rho,        
                    DUMMY_CAT //This must be the last one
};

class eventCategory{
  public:

    eventCategory(const std::string & aName, std::vector<eventCategory*> & categoryRejester){
      myName = aName;

      myId = categoryRejester.size();
      categoryRejester.push_back(this);

      std::cout<<__func__<<" "
      <<categoryRejester.size()
      <<" name: "<<myName
      <<std::endl;

      myWCategory = 0;
      myQCDCategory = 0;
      if(aName.find("_W")==std::string::npos) myWCategory = new eventCategory(aName+"_W",categoryRejester);
      if(aName.find("_QCD")==std::string::npos) myQCDCategory = new eventCategory(aName+"_QCD",categoryRejester);
    };

    unsigned int id() const {return myId;};

    const eventCategory * wControl() const {
      if(myWCategory) return myWCategory;
      else return this;
    };
    const eventCategory * qcdControl() const {
      if(myQCDCategory) return myQCDCategory;
      else return this;
    };

    std::string name() const {return myName;}

  private:

    std::string myName;
    unsigned int myId;
    eventCategory *myWCategory;
    eventCategory *myQCDCategory;

};

//eventCategory test1("jet0");
//eventCategory test2("boosted");

/*
//Signal categories
eventCategory jet0;
eventCategory boosted;
eventCategory vbf;
//Tight to Loose ratio categories
eventCategory tight_ss_jet0;
eventCategory tight_ss_boosted;
eventCategory tight_ss_vbf;

eventCategory loose_ss_jet0;
eventCategory loose_ss_boosted;
eventCategory loose_ss_vbf;

//antiiso control region
eventCategory antiIso_jet0;
eventCategory antiIso_boosted;
eventCategory antiIso_vbf;

//cp categories
eventCategory mu_pi;
eventCategory mu_rho;

eventCategory pi_pi;
eventCategory pi_rho;
eventCategory rho_rho;
*/

enum sysEffects{NOMINAL, NOMINAL_SVFIT,
		    TESUp, TESDown,
		    JESUp, JESDown,
		    M2TUp, M2TDown,
		    E2TUp, E2TDown,
		    J2TUp, J2TDown,
		    ZPtUp, ZPtDown,
		    TTUp, TTDown,
		    QCDSFUp, QCDSFDown,
		    WSFUp, WSFDown,
        DUMMY_SYS,
        //syst effects not implemented in data should be put at the end.
        //histograms will be the same as for NOMINAL
        ggUp, ggDown,
        ZmumuUp, ZmumuDown
      };
}
#endif
