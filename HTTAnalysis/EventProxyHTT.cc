#include "EventProxyHTT.h"

EventProxyHTT::EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHTT::~EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyHTT::init(std::vector<std::string> const& iFileNames){

	treeName_ = "m2n/event";

	EventProxyBase::init(iFileNames);

	fChain->SetBranchStatus("*",1);
	fChain->SetMakeClass(0);
       
	fChain->SetBranchAddress("npv",&npv);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////