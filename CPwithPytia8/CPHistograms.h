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

  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);

   private:
  
  virtual void defineHistograms();

  void plotHistograms(const std::string & sysType);

  ///Plot histogram hName+sysType for h, A, and Z on single canvas.
  void plot_HAZ_Histograms(const std::string & hName,
			   const std::string & sysType);

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
