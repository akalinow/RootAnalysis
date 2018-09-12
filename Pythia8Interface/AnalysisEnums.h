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

class eventCategory {
public:

eventCategory(const std::string & aName, std::vector<const eventCategory*> & categoryRejester){
        myName = aName;

        myId = categoryRejester.size();
        categoryRejester.push_back(this);

        myWSF = 0;
        myQCDEstimate = 0;
        myQCDSFNum = 0;
        myQCDSFDenom = 0;
        if(aName.find("_W")==std::string::npos) myWSF = new eventCategory(aName+"_W",categoryRejester);
        if(aName.find("_QCD")==std::string::npos) {
                myQCDEstimate = new eventCategory(aName+"_QCD",categoryRejester);
                myQCDSFNum = new eventCategory(aName+"_QCDSFNum",categoryRejester);
                myQCDSFDenom = new eventCategory(aName+"_QCDSFDenom",categoryRejester);
        }
};

unsigned int id() const {
        return myId;
};

const eventCategory * wSF() const {
        if(myWSF) return myWSF;
        else return this;
};
const eventCategory * qcdEstimate() const {
        if(myQCDEstimate) return myQCDEstimate;
        else return this;
};
const eventCategory * qcdSFNumerator() const {
        if(myQCDSFNum) return myQCDSFNum;
        else return this;
};
const eventCategory * qcdSFDenominator() const {
        if(myQCDSFDenom) return myQCDSFDenom;
        else return this;
};

const std::string & name() const {
        return myName;
}

private:

std::string myName;
unsigned int myId;
eventCategory *myWSF;
eventCategory *myQCDEstimate;
eventCategory *myQCDSFNum;
eventCategory *myQCDSFDenom;


};

enum sysEffects {NOMINAL,
                 TESUp, TESDown,
                 JESUp, JESDown,
                 M2TUp, M2TDown,
                 E2TUp, E2TDown,
                 DUMMY_SYS,
///Place systematic effects not affecting the SV calculation after DUMMY_SYS
///all quantities for following syst effects are calculated onfly, no need to rerun
///the ntuple making step.
                 J2TUp, J2TDown,
                 ZPtUp, ZPtDown,
                 TTUp, TTDown,
                 QCDSFUp, QCDSFDown,
                 WSFUp, WSFDown,
                 ggUp, ggDown,
                 ZmumuUp, ZmumuDown};
}
#endif
