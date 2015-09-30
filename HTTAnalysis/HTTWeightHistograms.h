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

  HTTWeightHistograms(TFileDirectory *myDir);

  HTTWeightHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTWeightHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
