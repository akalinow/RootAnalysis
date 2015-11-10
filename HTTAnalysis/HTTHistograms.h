
#ifndef HTTHistograms_h
#define HTTHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include "AnalysisHistograms.h"

class THStack;


class HTTHistograms: public AnalysisHistograms {
   public:

  HTTHistograms(std::string fileName="Histos.root", int opt=0);

  HTTHistograms(TFileDirectory *myDir);

  HTTHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);

  float getLumi();

  float getSampleNormalisation(const std::string & sampleName);
  TH1* getQCDbackground(std::string, int);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot stacked histograms for each contributing process.
  //varName - name of variable to be plotted,
  //selType - selection type, i.e. baseline, background estimation, etc.
  THStack* plotStack(std::string varName, int selType);

  //Plot a single histogram. One has to provide the full
  //histogram name, e.g. including h1D prefix.
  void plotAnyHistogram(const std::string & hName);


};

#endif
