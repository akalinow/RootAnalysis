//! File defining an interface (just static std::map) to get PropertyEnum members by using their std::string names.
/*!
  \author Rafal Maselek
  \date May 2018
  
  See also PropertyEnum.h
*/

#include "PropertyEnum.h"
	
const std::map<std::string, PropertyEnum> PropertyEnumString::enumMap={
	{"PDGId", PropertyEnum::PDGId},
    {"charge", PropertyEnum::charge},
    {"decayMode", PropertyEnum::decayMode},
	{"discriminator", PropertyEnum::discriminator}, 
	{"muonID", PropertyEnum::muonID}, 
	{"typeOfMuon", PropertyEnum::typeOfMuon},
	{"byCombinedIsolationDeltaBetaCorrRaw3Hits", PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits},
	{"byIsolationMVArun2v1DBoldDMwLTraw", PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw},
	{"againstElectronMVA5category", PropertyEnum::againstElectronMVA5category},
	{"dxy", PropertyEnum::dxy},
    {"dz", PropertyEnum::dz},
    {"SIP", PropertyEnum::SIP},
	{"tauID", PropertyEnum::tauID}, 
	{"combreliso", PropertyEnum::combreliso}, 
	{"leadChargedParticlePt", PropertyEnum::leadChargedParticlePt},
	{"isGoodTriggerType", PropertyEnum::isGoodTriggerType}, 
	{"FilterFired", PropertyEnum::FilterFired},	 
	{"L3FilterFired", PropertyEnum::L3FilterFired},
	{"L3FilterFiredLast", PropertyEnum::L3FilterFiredLast},
	{"mc_match", PropertyEnum::mc_match}, 
	{"rawPt", PropertyEnum::rawPt}, 
	{"area", PropertyEnum::area}, 
	{"PUJetID", PropertyEnum::PUJetID}, 
	{"jecUnc", PropertyEnum::jecUnc},
	{"Flavour", PropertyEnum::Flavour}, 
	{"bDiscriminator", PropertyEnum::bDiscriminator}, 
	{"bCSVscore", PropertyEnum::bCSVscore}, 
	{"PFjetID", PropertyEnum::PFjetID}, 
	{"NONE", PropertyEnum::NONE}
};