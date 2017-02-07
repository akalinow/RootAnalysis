#ifndef RootAnalysis_Tools_H
#define RootAnalysis_Tools_H

#include <string>
#include "Tools.h"
#include "AnalysisEnums.h"

namespace HTTAnalysis {

std::string categoryName(unsigned int iCategory){
        if(iCategory==(int)jet0_low) return "jet0_low";
        else if(iCategory==(int)jet0_high) return "jet0_high";
        else if(iCategory==(int)jet1_low) return "jet1_low";
        else if(iCategory==(int)jet1_high) return "jet1_high";
        else if(iCategory==(int)vbf_low) return "vbf_low";
        else if(iCategory==(int)vbf_high) return "vbf_high";
        else if(iCategory==(int)W) return "W";
        else if(iCategory==(int)TT) return "TT";
        else if(iCategory==(int)jet0) return "0jet";
        else if(iCategory==(int)boosted) return "boosted";
        else if(iCategory==(int)vbf) return "vbf";
        else if(iCategory==(int)CP_Pi) return "CP_Pi";
        else if(iCategory==(int)CP_Rho) return "CP_Rho";
        else if(iCategory==(int)wjets_jet0) return "wjets_0jet";
        else if(iCategory==(int)wjets_boosted) return "wjets_boosted";
        else if(iCategory==(int)wjets_vbf) return "wjets_vbf";
        else if(iCategory==(int)antiiso_jet0) return "antiiso_0jet";
        else if(iCategory==(int)antiiso_boosted) return "antiiso_boosted";
        else if(iCategory==(int)antiiso_vbf) return "antiiso_vbf";
        return "Unknown";
}

std::string systEffectName(unsigned int iSystEffect){
        if(iSystEffect==(int)NOMINAL) return "";
        else if(iSystEffect==(int)NOMINAL_SVFIT) return "";
        else if(iSystEffect==(int)TESUp) return "_CMS_shape_t_mt_13TeVUp";
        else if(iSystEffect==(int)TESDown) return "_CMS_shape_t_mt_13TeVDown";
        else if(iSystEffect==(int)JESUp) return "_CMS_scale_j_13TeVUp";
        else if(iSystEffect==(int)JESDown) return "_CMS_scale_j_13TeVDown";
        else if(iSystEffect==(int)M2TUp) return "_CMS_htt_ZLShape_mt_13TeVUp";
        else if(iSystEffect==(int)M2TDown) return "_CMS_htt_ZLShape_mt_13TeVDown";
        else if(iSystEffect==(int)E2TUp) return "_CMS_htt_ZEShape_et_13TeVUp";
        else if(iSystEffect==(int)E2TDown) return "_CMS_htt_ZEShape_et_13TeVDown";
        else if(iSystEffect==(int)J2TUp) return "_CMS_htt_jetToTauFake_13TeVUp";
        else if(iSystEffect==(int)J2TDown) return "_CMS_htt_jetToTauFake_13TeVDown";
        else if(iSystEffect==(int)ZPtUp) return "_CMS_htt_dyShape_13TeVUp";
        else if(iSystEffect==(int)ZPtDown) return "_CMS_htt_dyShape_13TeVDown";
        else if(iSystEffect==(int)TTUp) return "_CMS_htt_ttbarShape_13TeVUp";
        else if(iSystEffect==(int)TTDown) return "_CMS_htt_ttbarShape_13TeVDown";
        else if(iSystEffect==(int)QCDSFUp) return "_QCDSFUncert_mt_CAT_13TeVUp";
        else if(iSystEffect==(int)QCDSFDown) return "_QCDSFUncert_mt_CAT_13TeVDown";
        else if(iSystEffect==(int)WSFUp) return "_WSFUncert_mt_CAT_13TeVUp";
        else if(iSystEffect==(int)WSFDown) return "_WSFUncert_mt_CAT_13TeVDown";
        return "_Unknown";
}

}
#endif
