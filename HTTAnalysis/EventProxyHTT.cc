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
	
	fChain->SetBranchAddress("wevent",&wevent);
	fChain->SetBranchAddress("wpair",&wpair);
	fChain->SetBranchAddress("wtau",&wtau);
	fChain->SetBranchAddress("wmu",&wmu);
	
	fChain->SetBranchStatus("*",1);

	fChain->SetBranchStatus("npu_",1);
	fChain->SetBranchStatus("sample_",1);
	fChain->SetBranchStatus("genevtweight_",1);

	//fChain->SetBranchStatus("wpair",1);
	//fChain->SetBranchStatus("wtau",1);
	//fChain->SetBranchStatus("wmu",1);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
