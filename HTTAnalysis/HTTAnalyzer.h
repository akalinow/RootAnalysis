#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>
#include <vector>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"
#include "RooWorkspace.h"

#include "Analyzer.h"

class HTTHistograms;

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;

class HTTAnalyzer: public Analyzer{

 public:

  ///Copy from DataFormats/TauReco/interface/PFTauDecayMode.h
  enum hadronicTauDecayModes
  {
    tauDecay1ChargedPion0PiZero,
    tauDecay1ChargedPion1PiZero,  // rho (770 MeV) mediated)
    tauDecay1ChargedPion2PiZero,  // a1  (1.2 GeV) mediated
    tauDecay1ChargedPion3PiZero,  // contaminated or unmerged photo
    tauDecay1ChargedPion4PiZero,  // contaminated or unmerged photo
    tauDecay2ChargedPion0PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion1PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion2PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion3PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion4PiZero,  // extra track or un-recod track
    tauDecay3ChargedPion0PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion1PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion2PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion3PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion4PiZero,  // a1  (1.2 GeV) mediated
    tauDecaysElectron,
    tauDecayMuon,
    tauDecayOther                 // catch-all
  };

  enum muTauCategory{jet0_low, jet0_high,
		     jet1_low, jet1_high,
		     vbf_low, vbf_high,
		     jet0, boosted, vbf,
		     CP_Pi, CP_Rho,
		     wjets_jet0,
		     wjets_boosted, wjets_vbf,
		     antiiso_jet0,
		     antiiso_boosted, antiiso_vbf,
		     W, TT,
		     DUMMY //This must be the last one
  };

  HTTAnalyzer(const std::string & aName);

  virtual ~HTTAnalyzer();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *);

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  static std::string categoryName(unsigned int iCategory){
    if(iCategory==(int)HTTAnalyzer::jet0_low) return "jet0_low";
    else if(iCategory==(int)HTTAnalyzer::jet0_high) return "jet0_high";
    else if(iCategory==(int)HTTAnalyzer::jet1_low) return "jet1_low";
    else if(iCategory==(int)HTTAnalyzer::jet1_high) return "jet1_high";
    else if(iCategory==(int)HTTAnalyzer::vbf_low) return "vbf_low";
    else if(iCategory==(int)HTTAnalyzer::vbf_high) return "vbf_high";
    else if(iCategory==(int)HTTAnalyzer::W) return "W";
    else if(iCategory==(int)HTTAnalyzer::TT) return "TT";
    else if(iCategory==(int)HTTAnalyzer::jet0) return "0jet";
    else if(iCategory==(int)HTTAnalyzer::boosted) return "boosted";
    else if(iCategory==(int)HTTAnalyzer::vbf) return "vbf";
    else if(iCategory==(int)HTTAnalyzer::CP_Pi) return "CP_Pi";
    else if(iCategory==(int)HTTAnalyzer::CP_Rho) return "CP_Rho";
    else if(iCategory==(int)HTTAnalyzer::wjets_jet0) return "wjets_0jet";
    else if(iCategory==(int)HTTAnalyzer::wjets_boosted) return "wjets_boosted";
    else if(iCategory==(int)HTTAnalyzer::wjets_vbf) return "wjets_vbf";
    else if(iCategory==(int)HTTAnalyzer::antiiso_jet0) return "antiiso_0jet";
    else if(iCategory==(int)HTTAnalyzer::antiiso_boosted) return "antiiso_boosted";
    else if(iCategory==(int)HTTAnalyzer::antiiso_vbf) return "antiiso_vbf";
    return "Unknown";
  }

  static std::string systEffectName(unsigned int iSystEffect){
   if(iSystEffect==(int)sysEffects::NOMINAL) return "";
   else if(iSystEffect==(int)sysEffects::NOMINAL_SVFIT) return "";
   else if(iSystEffect==(int)sysEffects::TESUp) return "_CMS_shape_t_mt_13TeVUp";
   else if(iSystEffect==(int)sysEffects::TESDown) return "_CMS_shape_t_mt_13TeVDown";
   else if(iSystEffect==(int)sysEffects::JESUp) return "_CMS_scale_j_13TeVUp";
   else if(iSystEffect==(int)sysEffects::JESDown) return "_CMS_scale_j_13TeVDown";
   else if(iSystEffect==(int)sysEffects::M2TUp) return "_CMS_htt_ZLShape_mt_13TeVUp";
   else if(iSystEffect==(int)sysEffects::M2TDown) return "_CMS_htt_ZLShape_mt_13TeVDown";
   else if(iSystEffect==(int)sysEffects::E2TUp) return "_CMS_htt_ZLShape_et_13TeVUp";
   else if(iSystEffect==(int)sysEffects::E2TDown) return "_CMS_htt_ZLShape_et_13TeVDown";
   else if(iSystEffect==(int)sysEffects::J2TUp) return "_CMS_htt_jetToTauFake_13TeVUp";
   else if(iSystEffect==(int)sysEffects::J2TDown) return "_CMS_htt_jetToTauFake_13TeVDown";
   else if(iSystEffect==(int)sysEffects::ZPtUp) return "_CMS_htt_dyShape_13TeVUp";
   else if(iSystEffect==(int)sysEffects::ZPtDown) return "_CMS_htt_dyShape_13TeVDown";
   else if(iSystEffect==(int)sysEffects::TTUp) return "_CMS_htt_ttbarShape_13TeVUp";
   else if(iSystEffect==(int)sysEffects::TTDown) return "_CMS_htt_ttbarShape_13TeVDown";
   else if(iSystEffect==(int)sysEffects::QCDSFUp) return "_QCDSFUncert_mt_CAT_13TeVUp";
   else if(iSystEffect==(int)sysEffects::QCDSFDown) return "_QCDSFUncert_mt_CAT_13TeVDown";
   else if(iSystEffect==(int)sysEffects::WSFUp) return "_WSFUncert_mt_CAT_13TeVUp";
   else if(iSystEffect==(int)sysEffects::WSFDown) return "_WSFUncert_mt_CAT_13TeVDown";
   return "_Unknown";
 }

  ///Check it the event passes given category selections.
  ///Selections common to all categories (mu pt, tau Id etc.)
  ///are checked outside this method.
  bool passCategory(const HTTAnalyzer::muTauCategory & aCategory,
		    const sysEffects::sysEffectsEnum & aSysEffect=sysEffects::NOMINAL);

  ///Check it tau decay modes (GEN and RECO) match selected (hardcoded)
  ///decay mode.
  std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it
  static std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it.
  ///Method used when sample coding in TTree is not present.
  ///In this case a ROOT file name is used to decode the sample type.
  static std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy);

  ///Return sample name for DY. Name encoded jet bin, and decay mode.
  static std::string getDYSampleName(const EventProxyHTT & myEventProxy);

  //Return name sample name suffix for different particles matched to reconstructed tau
  static std::string getMatchingName(const EventProxyHTT & myEventProxy);

  ///Return pileup reweighting weight.
  ///Weight is calculatedon fly using the ration of nPU
  ///histograms for data and analyased sample.
  float getPUWeight(const EventProxyHTT & myEventProxy);

  ///Return event weight for systematic effects
  ///implemented by a global event weight.
  float getSystWeight(const sysEffects::sysEffectsEnum & aSystEffect=sysEffects::NOMINAL);

  ///Fill pulls between generator and various reco vertices.
  bool fillVertices(const std::string & sysType);

  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const std::string & hNameSuffix, float eventWeight,
			 const sysEffects::sysEffectsEnum & aSysEffect=sysEffects::NOMINAL);


  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for reconstructed
  ///mu+tau_h mode
  void fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);

  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for
  ///generator level taus for all decay modes.
  void fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);

  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter,
					    bool sgn=true);

  ///Return string encoding di-tau decay mode.
  ///The event can belong to more than one category
  std::vector<std::string> getTauDecayName(int decModeMinus, int decModePlus);

  ///Check if the decMode points to single prong hadronic tau decay
  bool isOneProng(int decMode);

  ///Check if the decMode points to leptonic tau decay
  bool isLepton(int decMode);

  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
					    float deltaR);

  ///Get lepton corrections
  float getLeptonCorrection(float eta, float pt, hadronicTauDecayModes tauDecayMode);

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  void getPreselectionEff(const EventProxyHTT & myEventProxy);

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Convert RooRealVar functions to histograms
  void initializeCorrections();

  ///Histograms storage.
  HTTHistograms *myHistos_;

  ///ROOT file with PU histogram
  TFile *puDataFile_, *puMCFile_;

  ///ROOT file containing current TTree
  TFile *ntupleFile_;

  ///Histogram with event counts filled during preselection step.
  TH1F *hStatsFromFile;

  ///Histograms with lepton corrections
  TH2F *h2DMuonIdCorrections, *h2DMuonIsoCorrections, *h2DMuonTrgCorrections;
  TH3F *h3DTauCorrections;

  ///Vector of PU histograms for MC samples
  std::vector<TH1F*> hPUVec_;

  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;

  ///Reconstructed objects selected for given event.
  HTTEvent aEvent;
  HTTPair aPair;
  std::string sampleName;

  HTTParticle aTau, aMuon, aMET;
  HTTParticle aGenHadTau, aGenMuonTau;
  HTTParticle aJet1, aJet2;
  std::vector<HTTParticle> aSeparatedJets;
  int nJets30;
  int nJetsInGap30;
  std::vector<bool> categoryDecisions;

  //cut on nPCA
  const float nPCAMin_;

};

#endif
