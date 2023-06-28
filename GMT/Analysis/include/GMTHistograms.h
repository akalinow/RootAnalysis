#ifndef GMTHistograms_H
#define GMTHistograms_H

#include "AnalysisHistograms.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

class GMTHistograms: public AnalysisHistograms {
public:
  
  GMTHistograms(std::string fileName="Histos.root", int opt=0);
  
  GMTHistograms(TDirectory *myDir);

  GMTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~GMTHistograms(); 

  virtual void finalizeHistograms();

  virtual std::string getTemplateName(const std::string& name);

  static const int color[6];
  static const int ptCutsGmt[4];
  static const int ptCutsOMTF[4];
  static const int ptCutsOMTFHigh[4];
  static const unsigned int nPtBins;
  static const float ptBins[36];
    
private:

  virtual void defineHistograms();


  TH1* Integrate(TH1 * histoD);

  TEfficiency * DivideErr(TH1D * h1, TH1D * h2,
                   const char * name="DivideErr",
                   const char * optErr ="");

  void plotEffPanel(const std::string & sysType, bool doHigh=false);

  void plotEffVsEta(const std::string & sysType);

  void plotEffVsVar(const std::string & sysType,
		    const std::string & varName);

  void plotVar(const std::string & sysType,
	       const std::string & varName);
  
  void plotEffPanelGen(const std::string & sysType, bool doHigh=false);


  void plotGMTVsOther(int iPt, std::string sysType="BMTF");

  void plotSingleHistogram(std::string hName);

  TH2F* makeRateWeights(TH2 *hOrig);
  TH1* getRateHisto(std::string sysType = "Vx",
		    std::string type = "Tot");
  void plotRate(std::string type);
  void plotEffVsRate(int iPtCut);
  
  float getEfficiency(TH2F *h2D, float ptCut);

  TEfficiency *getEfficiencyHisto(const std::string & hName);
 
  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
