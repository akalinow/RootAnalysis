#ifndef OMTFHistograms_H
#define OMTFHistograms_H

#include "AnalysisHistograms.h"

class TEfficiency;

class OMTFHistograms: public AnalysisHistograms {
public:
  
  OMTFHistograms(std::string fileName="Histos.root", int opt=0);
  
  OMTFHistograms(TDirectory *myDir);

  OMTFHistograms(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~OMTFHistograms(); 

  virtual void finalizeHistograms();

  virtual std::string getTemplateName(const std::string& name);

  static const std::vector<std::string> algos;
  static const std::vector<double> color;
  static const std::vector<double> iPtCuts;
  static const std::vector<double> ptBins;

  static constexpr int lhcNumberOfBunches = 2760; //TEST 2760 2345;
  static constexpr int lhcRevolutionFreq = 11264;
  static constexpr double lhcRate = lhcNumberOfBunches * lhcRevolutionFreq;
  
private:

  virtual void defineHistograms();

  TH1* Integrate(TH1 * histoD);

  void plotEffPanel(const std::string & sysType, const std::string & varName);

  void plotEffVsEta(const std::string & sysType);

  void plotEffVsVar(const std::string & sysType,
		    const std::string & varName);

  void plotVar(const std::string & sysType,
	       const std::string & varName);

  void plotEffType1VsType2(int iPt, std::string sysType1="OMTF", std::string sysType2="BMTF");

  void plotSingleHistogram(std::string hName);

  double vxMuRate(double pt_GeV) const;
  TH2F* makeRateWeights(TH2 *hOrig, const std::string & flavour);
  TH1* getRateHisto(std::string sysType = "Vx", std::string type = "Tot");
  void plotRate(std::string type);
    
  TEfficiency * getEfficiency(const std::string & hName);
 
  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
