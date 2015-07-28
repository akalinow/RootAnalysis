
#ifndef EVENTPROXYCPNtuple_h
#define EVENTPROXYCPNtuple_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"

#include "TBranch.h"
#include "TLorentzVector.h"
#include "TVector3.h"

   class EventProxyCPNtuple: public EventProxyBase{

   public:

      EventProxyCPNtuple();
      virtual ~EventProxyCPNtuple();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

      // Declaration of leaf types
      TLorentzVector  *p4Sum;
      TLorentzVector  *metNu, *met;
      TLorentzVector  *piMinus, *piPlus;
      TLorentzVector  *tauMinus, *tauPlus;
      TLorentzVector  *visTauMinus, *visTauPlus;
      Int_t           bosonId, decModeMinus, decModePlus;
      TVector3        *thePV, *svMinus, *svPlus;
      TVector3        *nPiMinus, *nPiPlus;
      Float_t         phi, rho, phi2, phi3;
      Float_t         yMinus, yPlus, yMinus2, yPlus2;
      Float_t         yMinusLab, yPlusLab, yMinusLab2, yPlusLab2;
      
   private:

      // List of branches
      TBranch        *b_p4Sum, *b_metNu, *b_met;
      TBranch        *b_piMinus, *b_piPlus;
      TBranch        *b_tauMinus, *b_tauPlus;
      TBranch        *b_visTauMinus, *b_visTauPlus;
      TBranch        *b_bosonId, *b_decModeMinus, *b_decModePlus;
      TBranch        *b_thePV, *b_svMinus, *b_svPlus;
      TBranch        *b_nPiMinus, *b_nPiPlus; 
      TBranch        *b_phi, *b_rho, *b_phi2, *b_phi3;
      TBranch        *b_yMinus, *b_yPlus, *b_yMinus2, *b_yPlus2;
      TBranch        *b_yMinusLab, *b_yPlusLab, *b_yMinusLab2, *b_yPlusLab2;
      
};
#endif
