#ifndef RootAnalysis_AnalysisEnums_H
#define RootAnalysis_AnalysisEnums_H

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
                    CP_Pi, CP_Rho,
                    wjets_jet0, wjets_boosted, wjets_vbf,
                    wjets_qcd_jet0, wjets_qcd_boosted, wjets_qcd_vbf,
                    qcd_jet0, qcd_boosted, qcd_vbf,
                    pipi, pirho, rhorho,
                    W, TT,
                    DUMMY_CAT //This must be the last one
};

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
		    DUMMY_SYS};
}
#endif
