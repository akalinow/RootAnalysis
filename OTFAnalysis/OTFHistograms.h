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


class OTFHistograms: public AnalysisHistograms {
   public:

  OTFHistograms(std::string fileName="Histos.root", int opt=0);

  OTFHistograms(TFileDirectory *myDir);

  OTFHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~OTFHistograms();

  void finalizeHistograms(int nRuns, float weight=1.0);

   private:

  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;


  static const unsigned int nPtBins;
  static const float ptBins[33];

};

#endif
