#include <sstream>

#include "HTTSynchNTuple.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::HTTSynchNTuple(const std::string & aName):HTTSynchNTupleBase(aName){ tmpName = "h1DXSignal";}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTuple::clone() const{

  HTTSynchNTuple* clone = new HTTSynchNTuple(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2){
  //Specific implementation for the mu+tau decay channel

  //Leg1: muon
  iso_1 =  leg2.getProperty(PropertyEnum::combreliso);
  /*
    trigweight_2;
    idisoweight_2;
  */

  //Leg2: tau
  iso_2 = leg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
 //FIXME: following properites should to be added
 /*  
     againstElectronLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronLooseMVA6);				
     againstElectronMediumMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronMediumMVA6);
     againstElectronTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronTightMVA6);
     againstElectronVLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVLooseMVA6);
     againstElectronVTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVTightMVA6);
     againstMuonLoose3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);
     againstMuonTight3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);
  */
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
