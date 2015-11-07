#ifndef SummaryAnalyzer_H
#define SummaryAnalyzer_H

#include <string>
#include <bitset>
#include <vector>

#include "Analyzer.h"

#include "TBits.h"
#include "TList.h"

class TH1F;
class TH1D;
class TH2F;
class TTree;
class TBranch;

class TreeAnalyzer;

class SummaryAnalyzer: public  Analyzer{

 public:

  SummaryAnalyzer(const std::string & aName);

  virtual ~SummaryAnalyzer();
  
  virtual void initialize(TFileDirectory&,
			  pat::strbitset *aSelections);  

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual void finalize();

  virtual void addBranch(TTree *tree);

  virtual void addCutHistos(TList *aList);

  TTree* getTree() { return mySelectionsTree_;};

  TList* getHistoList() { return &cutCounterHistos_;};


 private:

  Analyzer* clone() const;

  void fillEffHisto(std::string type);

  ///Remove bins corresponding to selections from other flaour,
  ///e.g ."Calo" from hito for "PF" selections.
  void cleanHisto(TH2F *hCutCounter);

  ///For CutCounters collapse the 2D histo to 1D
  ///i.e. integrate out the X axis.
  void clean2DHisto(TH2F *hCutCounter);

  TH1D * Integrate(TH1D * histoD);

  TFileDirectory* myDir_;

  TH1F *hCutNames_;

  TTree *mySelectionsTree_;

  TBranch *branchWeight_, *bitsBranch_;

  ///Types of the selections, eg. "PF", "Calo"
  std::vector<std::string> selectionFlavours_;

  std::vector<float> cutCounters_;

  TList cutCounterHistos_;


  float eventWeight_;

  TBits selectionWord_;
 
};

#endif
