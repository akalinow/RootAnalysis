#include "EventProxyHTT.h"

#include "TSystem.h"

#include <iostream>

EventProxyHTT::EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHTT::~EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyHTT::init(std::vector<std::string> const& iFileNames){

	treeName_ = "m2n/newevent";

	EventProxyBase::init(iFileNames);

	fChain->SetMakeClass(0);
		
	wevent = 0;//IMPORTNANT!!
	wpair = 0;//IMPORTNANT!!
	wtau = 0;//IMPORTNANT!!
	wmu = 0;//IMPORTNANT!!
	wjet = 0;//IMPORTNANT!!
	
	fChain->SetBranchAddress("wevent",&wevent);
	fChain->SetBranchAddress("wpair",&wpair);
	fChain->SetBranchAddress("wtau",&wtau);
	fChain->SetBranchAddress("wmu",&wmu);
	fChain->SetBranchAddress("wjet",&wjet);
	
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
