#ifndef RootAnalysis_HTTSynchNTuple_H
#define RootAnalysis_HTTSynchNTuple_H

#include <string>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "HTTSynchNTupleBase.h"

class EventProxyHTT;
class HTTHistograms;

class HTTSynchNTuple: public HTTSynchNTupleBase{

 public:

  HTTSynchNTuple(const std::string & aName);
  
  Analyzer* clone() const;

  virtual void fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2);

 protected:


 private:

};

#endif
