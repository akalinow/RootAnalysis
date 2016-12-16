#ifndef HTTHistograms_h
#define HTTHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

#include "HTTEvent.h"
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
  TH1F* getQCDbackground(unsigned int iCategory, std::string varName,
			 std::pair<float,float> wselOSCorrection =  std::pair<float,float>(1,0),
			 std::pair<float,float> wselSSCorrection =  std::pair<float,float>(1,0),
       unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT);

  ///Calculate scaling factor for the WJets MC
  ///Scaling factor is estimated in high Mt region.
  ///Other backgrounds are subtracted, basing on MC
  ///QCD contribution is neglected.
  std::pair<float,float> getWNormalisation(unsigned int iCategory, std::string selName,
    unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT);

  ///Calculate QCD OS/SS ratiousing non isolated events.
  std::pair<float,float> getQCDOStoSS(unsigned int iCategory,
				      std::pair<float,float> wselOSCorrection =  std::pair<float,float>(1,0),
				      std::pair<float,float> wselSSCorrection =  std::pair<float,float>(1,0),
              unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT);

   private:

  std::pair<float,float> wselOSCorrection;
  std::pair<float,float> wselSSCorrection;

  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  //Plot stacked histograms for each contributing process.
  ///iCategory - selection category to be plotted.
  ///iSystEffect - systematic effect type to be plotted
  ///selName - secondary type of selection (OS/SS/mt) used for background estimation
  //varName - name of variable to be plotted,
  THStack* plotStack(unsigned int iCategory,
		     std::string varName,
		     std::string selName = "OS",
		     unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL);

  void plotnPCA(const std::string & type);

  void plotVerticesPulls(const std::string & hName);

  void plotProfiles(const std::string & hName,
		    const std::string & sysType);

  void plotPhiDecayPlanes(const std::string& name);

  void plot_HAZ_Histograms(const std::string & hName,
			   const std::string & sysType);

  void plotCPhistograms(unsigned int iCategory);

  ///Return histogram for sum of all DY decay modes, and jet bins
  TH1F *get1D_DYJet_Histogram(const std::string& name);

  ///Return histogram for sum of all jet bins
  TH1F *get1D_WJet_Histogram(const std::string& name);

  ///Return histogram for sum VV processes
  TH1F *get1D_VV_Histogram(const std::string& name, std::string tauMatchSuffix = "");

  ///Return histogram for sum TT processes
  TH1F *get1D_TT_Histogram(const std::string& name, std::string tauMatchSuffix = "");

  ///Return histogram for sum single top processes
  TH1F *get1D_ST_Histogram(const std::string& name);

  ///Return sum of DY histograms. Sum can run over
  ///decay modes, jet bins, or both.
  TH1F *get1D_TauMatchJetSum(const std::string& name, bool sumDecayModes, bool sumJetBins = false);

  //Return sum of histograms where you substitute pattern with all sample names.
  //If tauMatchSuffix is specified, you will have sum of histos that match the suffix
  TH1F* get1D_SumPattern_Histogram(const std::string& name, std::string pattern, std::vector<std::string> sampleNames, std::string tauMatchSuffix = "");

  ///Return histogram for sum of all W/Z jet bins
  ///The results is scaled to 1/LO_xsection.
  TH1F *get1D_VJetSum(const std::string& name);

  ///Return histogram for sum of the EWK + 2jets samples
  TH1F *get1D_EWK2JetsSum(const std::string& name);

  ///Return histogram from nJets sample normalised by
  ///preselection/number of analysed events
  TH1F *getNormalised_NJet_Histogram(const std::string& hName);

  //Plot a single histogram. One has to provide the full
  //histogram name, e.g. including h1D prefix.
  void plotSingleHistogram(std::string hName);

  float muTauDYScale, mumuDYScale;
  float ttScale;

};

#endif
