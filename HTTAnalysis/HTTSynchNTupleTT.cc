#include <sstream>

#include "HTTSynchNTupleTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTupleTT::HTTSynchNTupleTT(const std::string & aName):HTTSynchNTupleBase(aName){ tmpName = "h1DXSignal";}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTupleTT::clone() const{

  HTTSynchNTupleTT* clone = new HTTSynchNTupleTT(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleTT::fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2){
  //Specific implementation for the mu+tau decay channel

  //Leg1: leading tau
  iso_1 = leg1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
  byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = leg1.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
  againstElectronLooseMVA6_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<10) == 1<<10;				
  againstElectronMediumMVA6_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<11) == 1<<11;
  againstElectronTightMVA6_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<12) == 1<<12;
  againstElectronVLooseMVA6_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<9) == 1<<9;
  againstElectronVTightMVA6_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<13) == 1<<13;
  againstMuonLoose3_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<7) == 1<<7;
  againstMuonTight3_1 = (leg2.getProperty(PropertyEnum::tauID) & 1<<8) == 1<<8;
  /*
    chargedIsoPtSum_1;
    decayModeFindingOldDMs_1;
    neutralIsoPtSum_1;
    puCorrPtSum_1;
    trigweight_1;
    idisoweight_1;
  */

  //Leg2: trailing tau
  iso_2 = leg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
  againstElectronLooseMVA6_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<10) == 1<<10;				
  againstElectronMediumMVA6_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<11) == 1<<11;
  againstElectronTightMVA6_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<12) == 1<<12;
  againstElectronVLooseMVA6_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<9) == 1<<9;
  againstElectronVTightMVA6_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<13) == 1<<13;
  againstMuonLoose3_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<7) == 1<<7;
  againstMuonTight3_2 = (leg2.getProperty(PropertyEnum::tauID) & 1<<8) == 1<<8;
  /*
    chargedIsoPtSum_2;
    decayModeFindingOldDMs_2;
    neutralIsoPtSum_2;
    puCorrPtSum_2;
    trigweight_2;
    idisoweight_2;
  */

  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
