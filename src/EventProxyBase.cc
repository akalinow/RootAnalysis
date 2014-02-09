#include "EventProxyBase.h"


EventProxyBase::EventProxyBase(std::vector<std::string> const& iFileNames):
  fileNames_(), chain_(0), treeName_("outTreePtOrd"), eventIndex_(0){

  chain_ = boost::shared_ptr<TChain>(new TChain(treeName_.c_str()));

  for (auto it= iFileNames.begin(), itEnd = iFileNames.end();it!=itEnd; ++it) {
	  chain_->Add(it->c_str(),-1);
  }

  //chain_->Print();
  accumulatedSize_ = chain_->GetEntries();
  Int_t cachesize = 10000000; //10 MBytes
  chain_->SetCacheSize(cachesize);
  chain_->AddBranchToCache("*",kTRUE);
  //chain_->SetParallelUnzip(kTRUE);

  chain_->SetBranchStatus("*",0);

}

EventProxyBase::~EventProxyBase(){}

//
// member functions
//

EventProxyBase const&
EventProxyBase::operator++(){

	chain_->GetEntry(eventIndex_++);

   return *this;
}

/** Go to the very first Event*/

EventProxyBase const&
EventProxyBase::toBegin()
{

	chain_->GetEntry(0);
	eventIndex_ = 0;
    return *this;
}


bool
EventProxyBase::atEnd() const{

  return eventIndex_>=accumulatedSize_;
}
