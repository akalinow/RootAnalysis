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

  HTTHistograms(TDirectory *myDir);

  HTTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HTTHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  std::string getTemplateName(const std::string& name);
  
  float getLumi();

  ///Return sample normalisation based only on
  ///luminosity and cross section.
  ///MC to DATA scaling factors should be applied
  ///on top of this normalisation.
  float getSampleNormalisation(std::string sampleName);

  ///Estimate QCD background using the SS/OS method.
  TH1F* getQCDbackground(std::string varName, std::string selName,
			 std::pair<float,float> wselOSCorrection =  std::pair<float,float>(1,0),
			 std::pair<float,float> wselSSCorrection =  std::pair<float,float>(1,0));

  ///Calculate scaling factor for the WJets MC
  ///SCaling factor is estimated in high Mt region.
  ///Other backgrounds are subtracted, basing on MC
  ///QCD contribution is neglected.
  std::pair<float,float> getWNormalisation(std::string selName);

  ///Calculate QCD OS/SS ratiousing non isolated events.
  std::pair<float,float> getQCDOStoSS(std::string selName,
				      std::pair<float,float> wselOSCorrection =  std::pair<float,float>(1,0),
				      std::pair<float,float> wselSSCorrection =  std::pair<float,float>(1,0));


   private:

  std::pair<float,float> wselOSCorrection;
  std::pair<float,float> wselSSCorrection;
  
  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot stacked histograms for each contributing process.
  //varName - name of variable to be plotted,
  //selName - selection type name. For baseline use empty string
  THStack* plotStack(std::string varName, std::string selName);

  void plotnPCA(const std::string & type);

  void plotVerticesPulls(const std::string & hName);

  void plotProfiles(const std::string & hName,
		    const std::string & sysType);
  
  void plotPhiDecayPlanes(const std::string& name);

  void plot_HAZ_Histograms(const std::string & hName,
			   const std::string & sysType);

  void plotCPhistograms(int nRuns, float weight);

  ///Return histogram for sum of all DY decay modes.
  TH1F *get1D_DY_Histogram(const std::string& name);

  ///Return histogram for sum of all WJet HT bins
  TH1F *get1D_WJet_Histogram(const std::string& name);

  ///Return histogram fro nJets sample normalised by
  ///preselection/number of analysed events
  TH1F *getNormalised_NJet_Histogram(const std::string& hName,
				     const std::string& jetsName);

  //Plot a single histogram. One has to provide the full
  //histogram name, e.g. including h1D prefix.
  void plotSingleHistogram(std::string hName);
  
  void plotSingleProfile(std::string hName);

  float muTauDYScale, mumuDYScale;

};

#endif
