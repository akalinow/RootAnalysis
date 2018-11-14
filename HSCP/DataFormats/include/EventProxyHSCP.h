#ifndef EVENTPROXYHSCP_h
#define EVENTPROXYHSCP_h

#include <string>
#include <typeinfo>
#include <vector>

#include "EventProxyBase.h"
#include "HSCPEvent.h"

#include "TBranch.h"

   class EventProxyHSCP: public EventProxyBase{

   public:

      EventProxyHSCP();
      virtual ~EventProxyHSCP();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

      HSCPEvent event;
  
      ///Reset the data members.
      void clear();
      
   };
#endif
