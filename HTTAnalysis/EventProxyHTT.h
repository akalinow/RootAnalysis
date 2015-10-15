#ifndef EVENTPROXYHTT_h
#define EVENTPROXYHTT_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
//#include "HTTEvent.h"

#include "TBranch.h"

   class EventProxyHTT: public EventProxyBase{

   public:

      EventProxyHTT();
      virtual ~EventProxyHTT();

      void init(std::vector<std::string> const& iFileNames);

      float puWeight;
      float pfJetPt;
      float ptL2, etaL2;
      float ptL1, etaL1;
      int tightestHPSMVAWP;
      int isPFMuon, isTightMuon;
      int muFlag, vetoEvent;
      float diTauVisMass;
      float diTauNSVfitMass;
      float MtLeg1MVA;
      float combRelIsoLeg1DBetav2;
      float diTauCharge;
      int pairIndex;
      ULong64_t run;
      float HLTx;
      float HLTmatch;
      float numPV;
      float sampleWeight;
      int genDecay;
      
   private:


   };
#endif
