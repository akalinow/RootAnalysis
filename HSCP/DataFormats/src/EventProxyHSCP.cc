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

  fChain->SetBranchAddress("event_nr",&event.eventId);
  fChain->SetBranchAddress("run_nr",&event.runId);
  
  fChain->SetBranchAddress("pt_1",&event.candidates[0].pt);
  fChain->SetBranchAddress("eta_1",&event.candidates[0].eta);
  fChain->SetBranchAddress("phi_1",&event.candidates[0].phi);

  fChain->SetBranchAddress("pt_2",&event.candidates[1].pt);
  fChain->SetBranchAddress("eta_2",&event.candidates[1].eta);
  fChain->SetBranchAddress("phi_2",&event.candidates[1].phi);
  

  fChain->SetBranchStatus("*",1);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHSCP::clear(){

  event.clear();

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
