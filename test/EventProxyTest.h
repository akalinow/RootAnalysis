#ifndef EVENTPROXYTEST_h
#define EVENTPROXYTEST_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"

#include "TBranch.h"

class Point{

 public:

  double x;
  double y;
  
};
//////////////////////////////////////////////
//////////////////////////////////////////////
   class EventProxyTest: public EventProxyBase{

   public:

      EventProxyTest();
      virtual ~EventProxyTest();

      void init(std::vector<std::string> const& iFileNames);

      Point *myPoint;

      double x;
      
   };
#endif
