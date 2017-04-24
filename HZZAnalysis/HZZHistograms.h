#ifndef HZZHistograms_h
#define HZZHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  pon, 9 lis 2015, 08:55:38 CET
//
//
#include "AnalysisHistograms.h"

class THStack;


class HZZHistograms: public AnalysisHistograms {
   public:

  HZZHistograms(std::string fileName="Histos.root", int opt=0);

  HZZHistograms(TDirectory *myDir);

  HZZHistograms(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HZZHistograms();

  void finalizeHistograms();

  virtual std::string getTemplateName(const std::string& name);

   private:

  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot a single histogram.
  void plotAnyHistogram(const std::string & hName);

};

#endif
