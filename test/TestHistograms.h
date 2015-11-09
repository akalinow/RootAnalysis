#ifndef TestHistograms_h
#define TestHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  pon, 9 lis 2015, 08:55:38 CET
//
//
#include "AnalysisHistograms.h"

class THStack;


class TestHistograms: public AnalysisHistograms {
   public:

  TestHistograms(std::string fileName="Histos.root", int opt=0);

  TestHistograms(TFileDirectory *myDir);

  TestHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~TestHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot a single histogram.
  void plotAnyHistogram(const std::string & hName);

};

#endif
