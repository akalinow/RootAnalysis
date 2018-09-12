//! Class for extracting data for further ML processing in python.
/*!
  \author Rafal Maselek
  \date May 2018
  
  This class extracts four-momenta of legs and jets from the analysis, together with some global parameters
  and particles' properties. Data is written to a TTree in "Summary" branch of the main analysis file. It
  can be then exported to python and TensorFlow.
*/
#ifndef RootAnalysis_MLAnalyzer_H
#define RootAnalysis_MLAnalyzer_H

#include "Analyzer.h"
#include "AnalysisEnums.h"
#include "MLObjectMessenger.h"
#include "TLorentzVector.h"
#include "HTTEvent.h"
#include <set>

class MLAnalyzer : public Analyzer
{

 public:
  //! Constructor with the name.
  MLAnalyzer(const std::string & aName, const std::string & aDecayMode);
  //! Destructor
  virtual ~MLAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections){

    //default is to do whatever the analyzer does
    filterEvent_ = false;
};
  

  virtual Analyzer* clone() const;
  
  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger);

  void performAnalysis(const EventProxyBase& iEvent, ObjectMessenger *aMessenger);

  virtual void finalize(){;};

  //! Clear the member field collections.
  virtual void clear();

  //! Add branches to the TTree.
  virtual void addBranch(TTree *tree);

  virtual void addCutHistos(TList *aList){;};

  inline bool filter() const{ return filterEvent_;};

  //! Name of the config file, where properties of particles are set to be extracted from analysis.
  std::string cfgFileName;
  //! Decay mode
  std::string decayMode;

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  int execution_counter_; /*!< Counts calls of analyze() */
  unsigned legs_no_; /*!< Number of legs in current analysis. Has to be constant during the execution. */
  unsigned jets_no_; /*!< Number of jets in current analysis. Has to be constant during the execution. */
  TTree *MLTree_; /*!< Stores the address of output TTree */


  //should this Analyzer be able to filter events
  bool filterEvent_;

  //! containers for four-momenta
  std::vector<double> legs_p4_;
  std::vector<double> jets_p4_;
  //! storage for params loaded from file and assignement to concrete legs and jets
  std::map<std::string, std::vector<double>> params_; /*!< current values of params */
  std::map<std::string, std::set<unsigned>> params_legs_; /*!< which legs have which params */
  std::map<std::string, std::set<unsigned>> params_jets_; /*!< which jets have which params */

  //! global parameters of the analylis
  int nJets30_; /*!< How many interesting jets were present in an event */
  float visMass_; 
  float higgsPT_;
  float betaScore_; /*!< Beta score of B-jet */
  float higgsMassTrans_; /*!< Transverse mass of Higgs */

  void extractP4(const TLorentzVector& v, const std::string destination, const unsigned no); /*!< copy p4 of jets and legs into jets_p4_ and legs_p4_ */
  void prepareVariables(const std::vector<const HTTParticle*>* legs, const std::vector<const HTTParticle*>* jets); /*!< Prepare vars that will be used for saving to TTree */
  void globalsHTT(const MLObjectMessenger* mess, const std::vector<const HTTParticle*>* legs, const HTTAnalysis::sysEffects* aSystEffect); /*!< Calculate global parameters for HTTAnalysis */ 
  void parseCfg(const std::string & cfgFileName); /*!< Parse the config file */
  void fixSelections(); /*!< Fix the numbers of legs for extraction in case when parameters for all legs are to be extracted */
  void extractParticleProperties(const std::vector<const HTTParticle*>* legs,  const std::vector<const HTTParticle*>* jets); /*!< Extract values of particle properties to variables in collections */



};



#endif
