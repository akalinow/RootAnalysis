#include "EventProxyHTT.h"

#include "TSystem.h"

#include <iostream>

EventProxyHTT::EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHTT::~EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyBase* EventProxyHTT::clone() const{

  return new EventProxyHTT();

}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void EventProxyHTT::init(std::vector<std::string> const& iFileNames){

  treeName_ = "HTauTauTree";

  EventProxyBase::init(iFileNames);
  //fChain->SetMakeClass(0);

  event = 0;
  pairs = 0;
  jets = 0;
  genLeptons = 0;

  fChain->SetBranchAddress("HTTEvent.",&event);
  fChain->SetBranchAddress("HTTPairCollection",&pairs);
  fChain->SetBranchAddress("HTTJetCollection",&jets);
  fChain->SetBranchAddress("HTTGenLeptonCollection",&genLeptons);

  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("HTTEvent*",1);
  fChain->SetBranchStatus("HTTPairCollection*",1);
  fChain->SetBranchStatus("HTTJetCollection*",1);
  fChain->SetBranchStatus("HTTGenLeptonCollection*",1);
  //fChain->SetBranchStatus("HTTLeptonCollection*",1);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHTT::enableBranches(){

  fChain->SetBranchStatus("*",1);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHTT::disableBranches(){

  fChain->SetBranchStatus("*",0);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHTT::clear(){

  if(event) event->clear();
  if(pairs) pairs->clear();
  if(jets) jets->clear();
  if(genLeptons) genLeptons->clear();

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
