#include "EventProxyHTT.h"

#include <iostream>

EventProxyHTT::EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHTT::~EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyHTT::init(std::vector<std::string> const& iFileNames){

	treeName_ = "outTreePtOrd";

	EventProxyBase::init(iFileNames);

	fChain->SetMakeClass(0);
	
	run = 0;
	fChain->SetBranchAddress("run",&run);
	fChain->SetBranchAddress("puWeight",&puWeight);
	fChain->SetBranchAddress("pfJetPt",&pfJetPt);
	fChain->SetBranchAddress("ptL1",&ptL1);
	fChain->SetBranchAddress("ptL2",&ptL2);
	fChain->SetBranchAddress("etaL1",&etaL1);
	fChain->SetBranchAddress("etaL2",&etaL2);
	fChain->SetBranchAddress("tightestHPSMVAWP",&tightestHPSMVAWP);
	fChain->SetBranchAddress("isPFMuon",&isPFMuon);
	fChain->SetBranchAddress("isTightMuon",&isTightMuon);
	fChain->SetBranchAddress("muFlag",&muFlag);
	fChain->SetBranchAddress("vetoEvent",&vetoEvent);
	fChain->SetBranchAddress("muFlag",&muFlag);
	fChain->SetBranchAddress("diTauVisMass",&diTauVisMass);
	fChain->SetBranchAddress("diTauNSVfitMass",&diTauNSVfitMass);
	fChain->SetBranchAddress("MtLeg1MVA",&MtLeg1MVA);
	fChain->SetBranchAddress("diTauVisMass",&diTauVisMass);
	fChain->SetBranchAddress("combRelIsoLeg1DBetav2",&combRelIsoLeg1DBetav2);
	fChain->SetBranchAddress("diTauCharge",&diTauCharge);
	fChain->SetBranchAddress("pairIndex",&pairIndex);
	fChain->SetBranchAddress("HLTx",&HLTx);	
	fChain->SetBranchAddress("HLTmatch",&HLTmatch);	
	fChain->SetBranchAddress("numPV",&numPV);
	fChain->SetBranchAddress("sampleWeight",&sampleWeight);
	fChain->SetBranchAddress("genDecay",&genDecay);	

	fChain->SetBranchStatus("*",0);	
	fChain->SetBranchStatus("puWeight",1);
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("pfJetPt",1);
	fChain->SetBranchStatus("ptL1",1);
	fChain->SetBranchStatus("ptL2",1);
	fChain->SetBranchStatus("etaL1",1);
	fChain->SetBranchStatus("etaL2",1);
	fChain->SetBranchStatus("tightestHPSMVAWP",1);
	fChain->SetBranchStatus("isPFMuon",1);
	fChain->SetBranchStatus("isTightMuon",1);
	fChain->SetBranchStatus("muFlag",1);
	fChain->SetBranchStatus("vetoEvent",1);
	fChain->SetBranchStatus("muFlag",1);
	fChain->SetBranchStatus("diTauVisMass",1);
	fChain->SetBranchStatus("diTauNSVfitMass",1);	
	fChain->SetBranchStatus("MtLeg1MVA",1);
	fChain->SetBranchStatus("diTauVisMass",1);
	fChain->SetBranchStatus("combRelIsoLeg1DBetav2",1);
	fChain->SetBranchStatus("diTauCharge",1);
	fChain->SetBranchStatus("pairIndex",1);
	fChain->SetBranchStatus("HLTx",1);
	fChain->SetBranchStatus("HLTmatch",1);	
	fChain->SetBranchStatus("numPV",1);
	fChain->SetBranchStatus("sampleWeight",1);
	fChain->SetBranchStatus("genDecay",1);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
