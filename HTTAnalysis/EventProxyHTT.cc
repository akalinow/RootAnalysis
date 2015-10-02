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
	TObjArray* myList = fChain->GetListOfFiles();
	for(unsigned int i=0;i<myList->GetSize();++i){
	  TFile *file = (TFile*)myList->At(i);
	  if(file) fFriendChain->Add(file->GetTitle(),-1);
	}
	fChain->AddFriend(fFriendChain);
	fChain->SetMakeClass(0);
	
	///Add weight friend TTree
	TChain *fFriendChainWeights = new TChain("Summary/tree");	
	TFile *file = new TFile("RootAnalysis_Weights.root");
	TTree *treeWeights;
	if(file->IsOpen()){
	  treeWeights = (TTree*)file->Get("Summary/tree");
	  //fFriendChainWeights->Add("RootAnalysis_Weights.root",-1);
	}
	/////////////
	//fChain->AddFriend(treeWeights,"weights");
	fChain->AddFriend(treeWeights);
	
	//fChain->GetFriend("weights")->Print();

	TObjArray *aObjArray = fChain->GetListOfBranches();
	for(unsigned int i=0;i<aObjArray->GetSize();++i){
	  TBranch *branch = (TBranch*)aObjArray->At(i);
	  //if(branch) std::cout<<branch->GetName()<<std::endl;
	}

	puWeight = -1;
	
	fChain->SetBranchAddress("npv",&npv);
	fChain->SetBranchAddress("run",&run);	
	fChain->SetBranchAddress("svfit",&svfit);
	fChain->SetBranchAddress("PUWeight",&puWeight);
	
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("npv",1);
	fChain->SetBranchStatus("run",1);
	fChain->SetBranchStatus("svfit",1);
	fChain->SetBranchStatus("PUWeight",1);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
