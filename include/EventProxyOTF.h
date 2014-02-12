#ifndef EVENTPROXYOTF_h
#define EVENTPROXYOTF_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
#include "TBranch.h"

 static const Int_t kMaxl1ObjectsOtf = 1;
 static const Int_t kMaxl1ObjectsGmt = 4;


   class EventProxyOTF: public EventProxyBase{

   public:

      EventProxyOTF();
      virtual ~EventProxyOTF();

      void init(std::vector<std::string> const& iFileNames);

      // Declaration of leaf types
     //EventData       *Events;
       UInt_t          fUniqueID;
       UInt_t          fBits;
       Float_t         weight;
       Float_t         pt;
       Float_t         eta;
       Float_t         phi;
       Float_t         phiHit;
       Float_t         etaHit;
       Int_t           charge;
       Int_t           l1ObjectsOtf_;
       UInt_t          l1ObjectsOtf_fUniqueID[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       UInt_t          l1ObjectsOtf_fBits[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Float_t         l1ObjectsOtf_pt[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Float_t         l1ObjectsOtf_eta[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Float_t         l1ObjectsOtf_phi[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Int_t           l1ObjectsOtf_bx[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Int_t           l1ObjectsOtf_q[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Int_t           l1ObjectsOtf_charge[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Int_t           l1ObjectsOtf_type[kMaxl1ObjectsOtf];   //[l1ObjectsOtf_]
       Int_t           l1ObjectsGmt_;
       UInt_t          l1ObjectsGmt_fUniqueID[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       UInt_t          l1ObjectsGmt_fBits[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Float_t         l1ObjectsGmt_pt[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Float_t         l1ObjectsGmt_eta[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Float_t         l1ObjectsGmt_phi[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Int_t           l1ObjectsGmt_bx[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Int_t           l1ObjectsGmt_q[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Int_t           l1ObjectsGmt_charge[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]
       Int_t           l1ObjectsGmt_type[kMaxl1ObjectsGmt];   //[l1ObjectsGmt_]


   private:

     // List of branches
     TBranch        *b_Events_fUniqueID;   //!
     TBranch        *b_Events_fBits;   //!
     TBranch        *b_Events_weight;   //!
     TBranch        *b_Events_pt;   //!
     TBranch        *b_Events_eta;   //!
     TBranch        *b_Events_phi;   //!
     TBranch        *b_Events_phiHit;   //!
     TBranch        *b_Events_etaHit;   //!
     TBranch        *b_Events_charge;   //!
     TBranch        *b_Events_l1ObjectsOtf_;   //!
     TBranch        *b_l1ObjectsOtf_fUniqueID;   //!
     TBranch        *b_l1ObjectsOtf_fBits;   //!
     TBranch        *b_l1ObjectsOtf_pt;   //!
     TBranch        *b_l1ObjectsOtf_eta;   //!
     TBranch        *b_l1ObjectsOtf_phi;   //!
     TBranch        *b_l1ObjectsOtf_bx;   //!
     TBranch        *b_l1ObjectsOtf_q;   //!
     TBranch        *b_l1ObjectsOtf_charge;   //!
     TBranch        *b_l1ObjectsOtf_type;   //!
     TBranch        *b_Events_l1ObjectsGmt_;   //!
     TBranch        *b_l1ObjectsGmt_fUniqueID;   //!
     TBranch        *b_l1ObjectsGmt_fBits;   //!
     TBranch        *b_l1ObjectsGmt_pt;   //!
     TBranch        *b_l1ObjectsGmt_eta;   //!
     TBranch        *b_l1ObjectsGmt_phi;   //!
     TBranch        *b_l1ObjectsGmt_bx;   //!
     TBranch        *b_l1ObjectsGmt_q;   //!
     TBranch        *b_l1ObjectsGmt_charge;   //!
     TBranch        *b_l1ObjectsGmt_type;   //!

};
#endif
