#ifndef EVENTPROXY_h
#define EVENTPROXY_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"

   class EventProxyBase{

   public:

      EventProxyBase(std::vector<std::string> const& iFileNames);
      virtual ~EventProxyBase();

      EventProxyBase const& operator++();

      // Go to the very first Event.
      EventProxyBase const& toBegin();
      virtual bool atEnd() const;

      bool isValid() const;
      operator bool() const;

      TFile* getTFile() const { return chain_->GetFile();}

      //Event *getEvent() const { return Event(chain_);}

      Long64_t size() const{ return accumulatedSize_; }

   private:

      std::vector<std::string> fileNames_;
      boost::shared_ptr<TChain> chain_;
      std::string treeName_;

      Long64_t eventIndex_;
      Long64_t accumulatedSize_;

};
#endif