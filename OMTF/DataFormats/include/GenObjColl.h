#ifndef GenObjColl_H
#define GenObjColl_H

#include "TObject.h"
#include "GenObj.h"

#include <vector>
#include <ostream>

class GenObjColl : public TObject {
public:
  GenObjColl(const std::vector<GenObj> & coll = std::vector<GenObj>()) : theColl(coll) {}
  operator const std::vector<GenObj> & () const { return theColl;}
  operator bool() const {return !theColl.empty(); }
  const std::vector<GenObj> & data() const { return theColl;}
private:
  std::vector<GenObj> theColl;
  friend std::ostream & operator<< (std::ostream &out, const GenObjColl & coll) {
    int i=-1; 
    for (const auto & obj : coll.theColl) { 
      out <<++i<<"_MUON: "<<obj; 
      if (i!=static_cast<int>(coll.theColl.size()-1) ) out <<std::endl; 
    }
    return out;
  }
  ClassDef(GenObjColl,1);
};


#endif
