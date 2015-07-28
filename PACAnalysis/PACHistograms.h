#ifndef TauAnalysis_FWLiteTools_JetVetoHistograms_h
#define TauAnalysis_FWLiteTools_JetVetoHistograms_h

// Base class for histogram managing.
//
// Original Author:  Artur Kalinowski
//         Created:  Wed Jul 22 12:56:54 CEST 2009
// $Id: JetVetoHistograms.h,v 1.6 2010/02/10 08:02:54 akalinow Exp $
//
//
#include "AnalysisHistograms.h"


class PACHistograms: public AnalysisHistograms {
   public:

  PACHistograms(std::string fileName="Histos.root", int opt=0);

  PACHistograms(TFileDirectory *myDir);

  PACHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~PACHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

  virtual bool fill2DHistogram(const std::string &name, float val1, float val2, float weight=1.0);

  static const int color[6];
  static const int ptCutsGmt[4];
  static const int ptCutsPAC[4];
  static const int ptCutsPACHigh[4];
  static const unsigned int nPtBins;
  static const float ptBins[33];


   private:

  virtual void defineHistograms();


  TH1* Integrate(TH1 * histoD);

  TH1D * DivideErr(TH1D * h1, TH1D * h2,
                   const char * name="DivideErr",
                   const char * optErr ="");

  void plotEffPanel(const std::string & sysType);

  void plotEffVsEta(const std::string & sysType);

  void plotEffVsVar(const std::string & sysType,
		  	  	  	const std::string & varName);

  void plotPACVsGmt(int iPt, std::string sysType="Gmt");

  TH2F* makeRateWeights(TH2 *hOrig);
  TH1* getRateHisto(std::string sysType = "Vx",
		    std::string type = "Tot");
  void plotRate(std::string type);
  void plotEffVsRate(int iPtCut);
  float getEfficiency(TH2F *h2D, float ptCut);

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
