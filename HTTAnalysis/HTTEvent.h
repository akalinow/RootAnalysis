#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_NEW_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_NEW_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include <vector>

//
// class declaration
//

class Wevent{

    int run_ = 0;
    float lumi_ = 0;
    long int event_ = 0;
    int nup_ = -999, npu_ = -999, npv_ = -999; 
    int paircount_ = -1;
    float genevtweight_ = 1;

  public:
    Wevent();
    ~Wevent();
    void run(int x){run_ = x;}
    void lumi(float x){lumi_ = x;}
    void event(long int x){event_ = x;}
    void nup(int x){nup_ = x;}
    void npu(int x){npu_ = x;}
    void npv(int x){npv_ = x;}
    void paircount(int x){paircount_ = x;}
    void genevtweight(float x){genevtweight_ = x;}

    int run(){return run_;}
    float lumi(){return lumi_;}
    long int event(){return event_;}
    int nup(){return nup_;}
    int npu(){return npu_;}
    int npv(){return npv_;}
    int paircount(){return paircount_;}
    float genevtweight(){return genevtweight_;}
};

enum tauidenum {
     decayModeFinding, 
     decayModeFindingNewDMs, 
     byCombinedIsolationDeltaBetaCorrRaw3Hits,
     byLooseCombinedIsolationDeltaBetaCorr3Hits,
     byMediumCombinedIsolationDeltaBetaCorr3Hits,
     byTightCombinedIsolationDeltaBetaCorr3Hits,
     chargedIsoPtSum,
     neutralIsoPtSum,
     puCorrPtSum,
     againstMuonLoose3,
     againstMuonTight3,
     againstElectronVLooseMVA5,
     againstElectronLooseMVA5,
     againstElectronMediumMVA5,
     againstElectronTightMVA5,
     againstElectronVTightMVA5,
     byIsolationMVA3newDMwoLTraw,
     byIsolationMVA3oldDMwoLTraw,
     byIsolationMVA3newDMwLTraw,
     byIsolationMVA3oldDMwLTraw,

};
class Wtau{
    
    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;
    float mt_ = -999;
    
    std::vector<float> tauID_ = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  public:
    Wtau();
    ~Wtau();
    void clear();

    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}
    void mt(float x){mt_ = x;}
    void tauID(tauidenum y, float x){tauID_[y] = x;}


    float pt(){return pt_;}
    float phi(){return phi_;}
    float eta(){return eta_;}
    float mass(){return mass_;}
    float charge(){return charge_;}
    float mt(){return mt_;}

    float tauID(tauidenum y){return tauID_[y];}

};
typedef std::vector<Wtau> WtauCollection;

class Wpair{
    
     float svfit_ = -999;
     float diq_ = -999;
     float pth_ = -999;
     float ptvis_ = -999;
     float m_vis_ = -999;
   
  public:
     Wpair();
     ~Wpair();
     
     void svfit(float x){svfit_ = x;}
     void diq(float x){diq_ = x;}
     void pth(float x){pth_ = x;}
     void ptvis(float x){ptvis_ = x;}
     void m_vis(float x){m_vis_ = x;}

     float svfit(){return svfit_;}
     float diq(){return diq_;}
     float pth(){return pth_;}
     float ptvis(){return ptvis_;}
     float m_vis(){return m_vis_;}

};
typedef std::vector<Wpair> WpairCollection;

enum triggersenum {
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1,
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1,
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2,
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2,
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3,
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3,
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1,
    HLT_IsoMu24_eta2p1_v1,
    HLT_IsoMu27_v1,
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2,
    HLT_IsoMu24_eta2p1_v2
};
class Wtriggers{

    std::vector<bool> triggers_ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
  public:
    Wtriggers();
    ~Wtriggers();

    void trigger(triggersenum y, bool x){triggers_[y] = x;}
    bool trigger(triggersenum y){return triggers_[y];}

};







#endif

