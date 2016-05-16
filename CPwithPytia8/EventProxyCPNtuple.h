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
      
      HTTEvent       *event;
      
   private:

      // List of branches
      TBranch        *b_p4Sum, *b_metNu, *b_met;
      TBranch        *b_piMinus, *b_piPlus;
      TBranch        *b_tauMinus, *b_tauPlus;
      TBranch        *b_visTauMinus, *b_visTauPlus;
      TBranch        *b_toyPiMinus, *b_toyPiPlus;
      TBranch        *b_toyNeutralMinus, *b_toyNeutralPlus;
      TBranch        *b_toyPiZeroMinus, *b_toyPiZeroPlus;
      TBranch        *b_toyTauMinus, *b_toyTauPlus;
      TBranch        *b_bosonId, *b_decModeMinus, *b_decModePlus;
      TBranch        *b_toyDecModeMinus, *b_toyDecModePlus;
      TBranch        *b_toyNChargedMinus, *b_toyNNeutralMinus, *b_toyNPiZeroMinus, 
	             *b_toyNChargedPlus,  *b_toyNNeutralPlus,  *b_toyNPiZeroPlus;
      TBranch        *b_thePV, *b_svMinus, *b_svPlus;
      TBranch        *b_nPiMinus, *b_nPiPlus; 
      TBranch        *b_toySvMinus, *b_toySvPlus, *b_toyNPiMinus, *b_toyNPiPlus; 
      TBranch        *b_phi, *b_rho, *b_phi2, *b_phiRho;
      TBranch        *b_yMinus, *b_yPlus, *b_yMinus2, *b_yPlus2;
      TBranch        *b_yMinusLab, *b_yPlusLab, *b_yMinusLab2, *b_yPlusLab2;
      TBranch        *b_yToyMinus, *b_yToyPlus;
      TBranch        *b_isoMinus, *b_isoPlus, *b_outerMinus, *b_outerPlus;
      
   };
#endif
