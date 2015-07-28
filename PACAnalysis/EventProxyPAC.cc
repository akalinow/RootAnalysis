#include "EventProxyPAC.h"


EventProxyPAC::EventProxyPAC(){}

EventProxyPAC::~EventProxyPAC(){}



void EventProxyPAC::init(std::vector<std::string> const& iFileNames){

	treeName_ = "efficiencyTree";

	EventProxyBase::init(iFileNames);

	events = 0;

	fChain->SetBranchStatus("*",1);
	fChain->SetMakeClass(0);
	fChain->SetBranchAddress("Events",&events);

}
