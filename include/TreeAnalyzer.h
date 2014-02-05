#ifndef  RootAnalysis_TreeAnalyzer_h
#define  RootAnalysis_TreeAnalyzer_h

#include <string>
#include <iostream>
#include <ostream>
#include <bitset>

#include "boost/unordered_map.hpp"

#include "TCollection.h"

#include "Analyzer.h"
#include "ObjectMessenger.h"

//#include "PhysicsTools/FWLite/interface/TFileService.h"

class TH1D;
class TProofOutputFile;
class FWLiteSummaryAnalyzer;
class ObjectMessenger;

class TreeAnalyzer: public Analyzer{

  friend class Analyzer;

public :

  ///Default conctructor takes the configuration cfg.py file name
  ///optional: pointer to the PROOF outputfile
  TreeAnalyzer(const std::string & name, 
		     const std::string & cfgFileName, 
		     TProofOutputFile * proofFile=0);

  ///Default destructor
  ~TreeAnalyzer();

  ///Initialisation method. The input cfg.py file is parsed here.
  void init(std::vector<Analyzer*> myAnalyzers);

  ///Finialisation of the analysis. Final statistics histograms is filled here.
  void finalize();

  ///Loop method, starting loop over the events.
  int loop();

  //Process the run info for the file
  virtual void processRunInfo(const edm::RunBase& aRun);

  ///Implementation of the abstract method
  virtual bool analyze(const edm::EventBase& iEvent) {return true;};

  virtual bool analyze(const edm::EventBase& iEvent,
		       const edm::RunBase& aRun);

  ///Method return the event weight for analyzed dataset.
  float getEventWeight() const {return eventWeight_;};

  ///Method returning the name of the configuration cfg.py file.
  std::string getCfgFileName() const {return cfgFileName_;};

  ///Method returning the name of the datasample analyzed.
  std::string getSampleName() const {return sampleName_;};

  ///Method returning the ROOT file whre the analysis histograms are stored. 
  TFile *getHistoFile() const { return &store_->file();};

  ///Method returning the ROOT file whre the analysis histograms are stored. 
  fwlite::TFileService *getTFileService() const { return store_;};

  ///Method to acces the selections word
  const pat::strbitset & getSelections() const {return *myStrSelections_;};

  ///Method normalising the histograms at the end of processing.
  ///Now only the histograms in Summary directory are normalised.
  void scaleHistograms();

 private:

   ///Method parsing the inpout cfg.py file.
   void parseCfg(const std::string & cfgFileName);

   ///Method reseting the state of analyzer for the new event.
   void clear();

   ///Parameter set for the analysis
   boost::shared_ptr<edm::ParameterSet> parameterSet_;

   ///Vector with list of analyzers.
   std::vector<Analyzer*> myAnalyzers_;

   ///Vector with list of TFileDirectory
   std::vector<TFileDirectory> myDirectories_;

   ///List of events to be read
   boost::unordered_map<std::string, bool> eventsToProcessHash_;


   ///Summary analyzer is a special one
   FWLiteSummaryAnalyzer *mySummary_;

   ///Path to the output file.
   std::string filePath_;

   ///Pointer to the ROOT file holding the histograms.
   fwlite::TFileService *store_;

   ///Histogram counting the number of event analyzed.
   ///Necessary for event counting when running on PROOF.
   TH1D *hStats_;

   ///List of the input EDM file names. This list is fetched from the
   ///cfg.py file.
   std::vector<std::string> fileNames_;

   ///Name of the analyzed sample.
   std::string sampleName_;

   ///Name of the input cfg.py file.
   std::string cfgFileName_;

   ///Number of events to be analyzed.
   int nEventsToAnalyze_;   

   ///Number of initial and final events in the
   ///processing of RECO/AOD with possible filtering, that
   ///produced the EDM filed read with this code
   unsigned int initialRecoEvents_, finalRecoEvents_;

   ///current run number
   unsigned int currentRun_;

   ///analyzed runs counter
   unsigned int nRunsAnalyzed_;

   ///Map of the selections. Selection name is the key in the map.
   pat::strbitset *myStrSelections_;


   ///Object holding information passed between analyzers
   ObjectMessenger* myObjMessenger_;

};

#endif


