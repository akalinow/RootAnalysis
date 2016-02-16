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

  treeName_ = "m2n/eventTree";
  
  EventProxyBase::init(iFileNames);
  fChain->SetMakeClass(0);
  
  wevent = 0;//IMPORTNANT!!
  wpair = 0;//IMPORTNANT!!
  wtau = 0;//IMPORTNANT!!
  wtauGen = 0;//IMPORTNANT!!
  wmu = 0;//IMPORTNANT!!
  wjet = 0;//IMPORTNANT!!
  wmet = 0;//IMPORTNANT!!
  
  fChain->SetBranchAddress("wevent",&wevent);
  fChain->SetBranchAddress("wpair",&wpair);
  fChain->SetBranchAddress("wtau",&wtau);
  fChain->SetBranchAddress("wtauGen",&wtauGen);
  fChain->SetBranchAddress("wmu",&wmu);
  fChain->SetBranchAddress("wjet",&wjet);
  fChain->SetBranchAddress("wmet",&wmet);
  
  fChain->SetBranchStatus("*",1);
  
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
  fChain->SetBranchStatus("wpair",1);  
  fChain->SetBranchStatus("sample_",1);
  fChain->SetBranchStatus("genevtweight_",1);
  fChain->SetBranchStatus("npu_",1);
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void  EventProxyHTT::clear(){
  
  wpair->clear();
  wtau->clear();
  wmu->clear();
  wjet->clear();
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
