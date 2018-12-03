#include "EventProxyHSCP.h"

#include "TSystem.h"

#include <iostream>

EventProxyHSCP::EventProxyHSCP(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHSCP::~EventProxyHSCP(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyBase* EventProxyHSCP::clone() const{

  return new EventProxyHSCP();

}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void EventProxyHSCP::init(std::vector<std::string> const& iFileNames){

  treeName_ = "demo/HscpTree";

  EventProxyBase::init(iFileNames);

  fChain->SetBranchAddress("event_num",&event.eventId);
  fChain->SetBranchAddress("run_num",&event.runId);
  
  fChain->SetBranchAddress("NHSCPs",&event.numberOfCandidates);
  fChain->SetBranchAddress("hscp_pt_trkref",&event.pt);
  fChain->SetBranchAddress("hscp_track_eta",&event.eta);
  fChain->SetBranchAddress("hscp_track_phi",&event.phi);

  fChain->SetBranchStatus("*",1);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHSCP::clear(){

  event.clear();

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
