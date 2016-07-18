#ifndef EVENTPROXYHTT_h
#define EVENTPROXYHTT_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
#include "HTTEvent.h"

#include "TBranch.h"

   class EventProxyHTT: public EventProxyBase{

   public:

      EventProxyHTT();
      virtual ~EventProxyHTT();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

      HTTEvent *event;
      std::vector<HTTPair>  *pairs;
      std::vector<HTTParticle>  *jets;
      std::vector<HTTParticle>  *genLeptons;
  
      ///Enable branches to be read
      void enableBranches();

      ///Disable selected branches 
      void disableBranches();

      ///Reset the data members.
      void clear();
      
   };
#endif
