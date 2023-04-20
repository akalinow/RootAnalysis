#ifndef RecoMuon_H
#define RecoMuon_H 

#include "TObject.h"
#include "RecoMuonObj.h"

#include <vector>
#include <ostream>

class RecoMuon : public TObject {
public:
	RecoMuon(const std::vector<RecoMuonObj> & coll = std::vector<RecoMuonObj>()) : theColl(coll) {}
	operator const std::vector<RecoMuonObj> & () const { return theColl;}
	operator bool() const {return !theColl.empty(); }
        const std::vector<RecoMuonObj> & data() const { return theColl;}
private:
	std::vector<RecoMuonObj> theColl;
        friend std::ostream & operator<< (std::ostream &out, const RecoMuon & coll) {
             int i=-1;
             for (const auto & obj : coll.theColl) {
                  out <<++i<<"_MUON: "<<obj;
                  if (i!=static_cast<int>(coll.theColl.size()-1) ) out <<std::endl;
              }
            return out;
         }
        ClassDef(RecoMuon,1);
};


#endif




