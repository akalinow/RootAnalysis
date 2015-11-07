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

      float puWeight;

      Wevent *wevent;
      std::vector<Wpair>  *wpair;
      std::vector<Wtau>  *wtau;
      std::vector<Wmu>  *wmu;
      
   };
#endif
