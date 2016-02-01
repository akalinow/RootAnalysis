#ifndef EVENTPROXYCPNtuple_h
#define EVENTPROXYCPNtuple_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
#include "HTTEvent.h"

#include "TBranch.h"

   class EventProxyCPNtuple: public EventProxyBase{

   public:

      EventProxyCPNtuple();
      virtual ~EventProxyCPNtuple();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

      // Declaration of leaf types
      HTTEvent       *event;
      
   private:

      // List of branches
      TBranch        *b_Events;

   };
#endif
