#include <iostream>

#include "EventProxyCPNtuple.h"

EventProxyCPNtuple::EventProxyCPNtuple(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyCPNtuple::~EventProxyCPNtuple(){}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
EventProxyBase* EventProxyCPNtuple::clone() const{

  return new EventProxyCPNtuple();
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyCPNtuple::init(std::vector<std::string> const& iFileNames){

	treeName_ = "genAna/hTTCPTree";

	EventProxyBase::init(iFileNames);

	std::cout<<"fChain: "<<fChain<<std::endl;

	fChain->SetBranchStatus("*",1);
	fChain->SetMakeClass(0);
	fChain->SetBranchAddress("p4Sum.", &p4Sum, &b_p4Sum);
	fChain->SetBranchAddress("metNu.", &metNu, &b_metNu);
	fChain->SetBranchAddress("met.", &met, &b_met);
	fChain->SetBranchAddress("piMinus.", &piMinus, &b_piMinus);
	fChain->SetBranchAddress("piPlus.", &piPlus, &b_piPlus);
	fChain->SetBranchAddress("tauMinus.", &tauMinus, &b_tauMinus);
	fChain->SetBranchAddress("tauPlus.", &tauPlus, &b_tauPlus);
	fChain->SetBranchAddress("visTauMinus.", &visTauMinus, &b_visTauMinus);
	fChain->SetBranchAddress("visTauPlus.", &visTauPlus, &b_visTauPlus);
	//
	fChain->SetBranchAddress("bosonId", &bosonId, &b_bosonId);
	fChain->SetBranchAddress("decModeMinus", &decModeMinus, &b_decModeMinus);
	fChain->SetBranchAddress("decModePlus", &decModePlus, &b_decModePlus);
	//
	fChain->SetBranchAddress("thePV.", &thePV, &b_thePV);
	fChain->SetBranchAddress("svMinus.", &svMinus, &b_svMinus);
	fChain->SetBranchAddress("svPlus.", &svPlus, &b_svPlus);
	fChain->SetBranchAddress("nPiMinus.", &nPiMinus, &b_nPiMinus);
	fChain->SetBranchAddress("nPiPlus.", &nPiPlus, &b_nPiPlus);
	//
	fChain->SetBranchAddress("phi", &phi, &b_phi);
	fChain->SetBranchAddress("rho", &rho, &b_rho);
	fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
	//fChain->SetBranchAddress("phi3", &phi3, &b_phi3);
	//
	fChain->SetBranchAddress("yMinus", &yMinus, &b_yMinus);
	fChain->SetBranchAddress("yPlus", &yPlus, &b_yPlus);
	fChain->SetBranchAddress("yMinus2", &yMinus2, &b_yMinus2);
	fChain->SetBranchAddress("yPlus2", &yPlus2, &b_yPlus2);
	fChain->SetBranchAddress("yMinusLab", &yMinusLab, &b_yMinusLab);
	fChain->SetBranchAddress("yPlusLab", &yPlusLab, &b_yPlusLab);
	fChain->SetBranchAddress("yMinusLab2", &yMinusLab2, &b_yMinusLab2);
	fChain->SetBranchAddress("yPlusLab2", &yPlusLab2, &b_yPlusLab2);


	std::cout<<"Here 0"<<std::endl;
	fChain->GetEntry(0);
	std::cout<<"Here 1"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
