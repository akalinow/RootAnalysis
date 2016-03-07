#ifndef HWvsEMULHistograms_h
#define HWvsEMULHistograms_h

//
// Original Author:  Artur Kalinowski
//         Created:  Wed Jul 22 12:56:54 CEST 2009
// $Id: JetVetoHistograms.h,v 1.6 2010/02/10 08:02:54 akalinow Exp $
//
//
#include "AnalysisHistograms.h"

class HWvsEMULHistograms: public AnalysisHistograms {
   public:
  void pieceHistogramsTogether();
  HWvsEMULHistograms(std::string fileName="Histos.root", int opt=0);

  HWvsEMULHistograms(TFileDirectory *myDir);

  HWvsEMULHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours);

  virtual ~HWvsEMULHistograms(); 

  void finalizeHistograms(int nRuns, float weight=1.0);
  
  virtual void finalizeHistograms();

  virtual bool fill1DHistogram(const std::string &name, float val1, float weight=1.0);

  virtual bool fill2DHistogram(const std::string &name, float val1, float val2, float weight=1.0);

  void plotPhi();
  void plotHWvsEMUL(std::string type);
  void plotVariable(std::string type);
  void plotDelta(std::string type);
  
  
   private:

  virtual void defineHistograms();

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

};

#endif
