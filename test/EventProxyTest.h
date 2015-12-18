#ifndef EVENTPROXYTEST_h
#define EVENTPROXYTEST_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"

#include "TBranch.h"

struct Point {
    Double_t x,y;
};
//////////////////////////////////////////////
//////////////////////////////////////////////
   class EventProxyTest: public EventProxyBase{

   public:

      EventProxyTest();
      virtual ~EventProxyTest();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

      Double_t x,y;
      
   };
#endif
