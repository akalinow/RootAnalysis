#ifndef HTTHistograms_h
#define HTTHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

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

  ///Return sample normalisation based only on
  ///luminosity and cross section.
  ///MC to DATA scaling factors should be applied
  ///on top of this normalisation.
  float getSampleNormalisation(std::string sampleName);

  ///Estimate QCD background using the SS/OS method.
  TH1F* getQCDbackground(std::string varName, std::string selName);

  ///Calculate scaling factor for the WJets MC
  ///SCaling factor is estimated in high Mt region.
  ///Other backgrounds are subtracted, basing on MC
  ///QCD contribution is neglected.
  std::pair<float,float> getWNormalisation(std::string selName);

  ///Calculate QCD OS/SS ratiousing non isolated events.
  std::pair<float,float> getQCDOStoSS(std::string selName);

   private:
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot stacked histograms for each contributing process.
  //varName - name of variable to be plotted,
  //selName - selection type name. For baseline use empty string
  THStack* plotStack(std::string varName, std::string selName);

  void plotPhiDecayPlanes(const std::string & sysType,
			  const std::string & cosType);

  ///Return histogram for sum of all DY decay modes.
  TH1F *get1D_DY_Histogram(const std::string& name);

  ///Return histogram for sum of all WJet HT bins
  TH1F *get1D_WJet_Histogram(const std::string& name);

  //Plot a single histogram. One has to provide the full
  //histogram name, e.g. including h1D prefix.
  void plotSingleHistogram(std::string hName);


};

#endif
