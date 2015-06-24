#ifndef CPHistograms_h
#define CPHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  śro, 24 cze 2015, 13:30:00 CEST
//
//
#include "AnalysisHistograms.h"


class CPHistograms: public AnalysisHistograms {
   public:

  CPHistograms(std::string fileName="Histos.root", int opt=0);

  CPHistograms(TFileDirectory *myDir);

  CPHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~CPHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill2DHistogram(const std::string &name, float val1, float val2, float weight=1.0);

   private:

  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
