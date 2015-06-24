
#ifndef EVENTPROXYCP_h
#define EVENTPROXYCP_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
//#include "EventData.h"

#include "TBranch.h"

   class EventProxyCP: public EventProxyBase{

   public:

      EventProxyCP();
      virtual ~EventProxyCP();

      void init(std::vector<std::string> const& iFileNames);

      // Declaration of leaf types
      //EventData       *events;


   private:

     // List of branches
     TBranch        *b_Events;

};
#endif
