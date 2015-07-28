#ifndef EVENTPROXYOTF_h
#define EVENTPROXYOTF_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
#include "EventData.h"

#include "TBranch.h"


   class EventProxyOTF: public EventProxyBase{

   public:

      EventProxyOTF();
      virtual ~EventProxyOTF();

      void init(std::vector<std::string> const& iFileNames);
      
      // Declaration of leaf types
      EventData       *events[128]{};

};
#endif
