#ifndef HTTWeightHistograms_h
#define HTTWeightHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include "AnalysisHistograms.h"


class HTTWeightHistograms: public AnalysisHistograms {
   public:

  HTTWeightHistograms(std::string fileName="Histos.root", int opt=0);

  HTTWeightHistograms(TDirectory *myDir);

  HTTWeightHistograms(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTWeightHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  std::string getTemplateName(const std::string& name);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
