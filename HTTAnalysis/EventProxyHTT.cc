#include "EventProxyHTT.h"

#include <iostream>

EventProxyHTT::EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyHTT::~EventProxyHTT(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyHTT::init(std::vector<std::string> const& iFileNames){

	treeName_ = "m2n/event";

	EventProxyBase::init(iFileNames);

	TChain *fFriendChain = new TChain("m2n/pair");
	TChain *fFriendChain1 = new TChain("m2n/weight");
	TObjArray* myList = fChain->GetListOfFiles();
	for(unsigned int i=0;i<myList->GetSize();++i){
	  TFile *file = (TFile*)myList->At(i);
	  if(file){
	    fFriendChain->Add(file->GetTitle(),-1);
	    fFriendChain1->Add(file->GetTitle(),-1);	
	  }
	}
	fChain->AddFriend(fFriendChain);
	fChain->AddFriend(fFriendChain1);
	fChain->SetMakeClass(0);
	
	///Add weight friend TTree
	TChain *fFriendChainWeights = new TChain("Summary/tree");	
	TFile *file = new TFile("RootAnalysis_Weights.root");
	TTree *treeWeights;
	if(file->IsOpen() && file->FindObjectAny("tree")){
	  treeWeights = (TTree*)file->Get("Summary/tree");
	  fChain->AddFriend(treeWeights);
	}
	//////////	
	puWeight = -1;
	genWeight = 1.0;
	
	fChain->SetBranchAddress("npv",&npv);
	fChain->SetBranchAddress("run",&run);	
	fChain->SetBranchAddress("svfit",&svfit);
	fChain->SetBranchAddress("PUWeight",&puWeight);
	fChain->SetBranchAddress("generator",&genWeight);

	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("npv",1);
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("svfit",1);
	fChain->SetBranchStatus("PUWeight",1);
	fChain->SetBranchStatus("generator",1);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
