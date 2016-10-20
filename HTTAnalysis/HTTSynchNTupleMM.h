#ifndef RootAnalysis_HTTSynchNTupleMM_H
#define RootAnalysis_HTTSynchNTupleMM_H

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

class HTTSynchNTupleMM: public HTTSynchNTupleBase{

 public:

  HTTSynchNTupleMM(const std::string & aName);
  
  Analyzer* clone() const;

  void fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2);

 protected:


 private:

};

#endif
