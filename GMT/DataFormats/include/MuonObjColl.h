#ifndef MuonObjColl_H
#define MuonObjColl_H 

#include "TObject.h"
#include "MuonObj.h"

#include <vector>
#include <ostream>

class MuonObjColl : public TObject {
public:
	MuonObjColl(const std::vector<MuonObj> & coll = std::vector<MuonObj>()) : theColl(coll) {}
	operator const std::vector<MuonObj> & () const { return theColl;}
	operator bool() const {return !theColl.empty(); }
        const std::vector<MuonObj> & getMuonObjs() const { return theColl; }
        const std::vector<MuonObj> & data() const { return theColl;}
private:
	std::vector<MuonObj> theColl;
        friend std::ostream & operator<< (std::ostream &out, const MuonObjColl & coll) {
             int i=-1;
             for (const auto & obj : coll.theColl) {
                  out <<++i<<"_MUON: "<<obj;
                  if (i!=static_cast<int>(coll.theColl.size()-1) ) out <<std::endl;
              }
            return out;
         }
        ClassDef(MuonObjColl,1);
};


#endif




