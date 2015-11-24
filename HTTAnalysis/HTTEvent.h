#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent1_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent1_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include <vector>

class DiTauData{
public:
///Data common for generator and reconstruction levels.
int decModeMinus_, decModePlus_;
TVector3 thePV_;
TVector3 svMinus_, svPlus_;
TVector3 nPiPlus_, nPiMinus_;
TVector3 nPiPlusAODvx_, nPiPlusGenvx_, nPiPlusRefitvx_;
TVector3 nPiMinusAODvx_, nPiMinusGenvx_, nPiMinusRefitvx_;
TLorentzVector p4Sum_;
TLorentzVector piMinus_, piPlus_;
TLorentzVector tauMinus_, tauPlus_;
TLorentzVector visTauMinus_, visTauPlus_;
///Reconstruction level only data
TVector3 pfPV_, pt2PV_, refitPfPV_, refitPfPVNoBS_;
bool isRefit_;
int nTracksInRefit_;
int pfPVIndex_, pt2PVindex_;
///Generator level only data
///Reset data members to default values.
void clear();
};
class HTTEvent {
public:
HTTEvent();
~HTTEvent();
///Reset data members to default values.
void clear();
float run_, lumi_, event_;
int bosonId_;
DiTauData genEvent_, recoEvent_;
};


class Wevent{

    int run_ = 0;
    float lumi_ = 0;
    long int event_ = 0;
    int nup_ = -999, npu_ = -999, npv_ = -999; 
    int paircount_ = -1;
    float genevtweight_ = 1;
    int sample_ = -1;

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
    void sample(int x){sample_ = x;}

    int run(){return run_;}
    float lumi(){return lumi_;}
    long int event(){return event_;}
    int nup(){return nup_;}
    int npu(){return npu_;}
    int npv(){return npv_;}
    int paircount(){return paircount_;}
    float genevtweight(){return genevtweight_;}
    int sample(){return sample_;}
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
    
    std::vector<bool> tauID_;//  = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1, -1 ,-1, -1, -1};

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


class Wmu{
    
    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;
    float mt_ = -999;
    float d0_ = -999;
    float dz_ = -999;
    float isLooseMuon_ = -999;
    float isTightMuon_ = -999;
    float isHighPtMuon_ = -999;
    float isMediumMuon_ = -999;
    float isTightnovtxMuon_ = -999;
    float iso_ = -999;

  public:
    Wmu();
    ~Wmu();
    void clear();

    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}
    void mt(float x){mt_ = x;}
    void d0(float x){d0_ = x;}
    void dz(float x){dz_ = x;}
    void isLooseMuon(float x){isLooseMuon_ = x;}
    void isTightMuon(float x){isTightMuon_ = x;}
    void isHighPtMuon(float x){isHighPtMuon_ = x;}
    void isMediumMuon(float x){isMediumMuon_ = x;}
    void isTightnovtxMuon(float x){isTightnovtxMuon_ = x;}
    void iso(float x){iso_ = x;}

    float pt(){return pt_;}
    float phi(){return phi_;}
    float eta(){return eta_;}
    float mass(){return mass_;}
    float charge(){return charge_;}
    float mt(){return mt_;}
    float d0(){return d0_;}
    float dz(){return dz_;}
    float isLooseMuon(){return isLooseMuon_;}
    float isTightMuon(){return isTightMuon_;}
    float isHighPtMuon(){return isHighPtMuon_;}
    float isMediumMuon(){return isMediumMuon_;}
    float isTightnovtxMuon(){return isTightnovtxMuon_;}
    float iso(){return iso_;}


};
typedef std::vector<Wmu> WmuCollection;

class Welectron{

    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;

  public:
    Welectron();
    ~Welectron();
    void clear();
    
    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}
 
    float pt(){return pt_;}
    float phi(){return phi_;}
    float eta(){return eta_;}
    float mass(){return mass_;}
    float charge(){return charge_;}
};
typedef std::vector<Welectron> WelectronCollection;

class Wmet{


    float metpx_ = -999;
    float metpt_ = -999;
    float metphi_ = -999;
    float mvacov00_ = -999;
    float mvacov01_ = -999;
    float mvacov10_ = -999;
    float mvacov11_ = -999;
 
  public:

    void metpx(float x){metpx_ = x;}
    void metpt(float x){metpt_ = x;}
    void metphi(float x){metphi_ = x;}
    void mvacov00(float x){mvacov00_ = x;}
    void mvacov01(float x){mvacov01_ = x;}
    void mvacov10(float x){mvacov10_ = x;}
    void mvacov11(float x){mvacov11_ = x;}

    float metpx(){return metpx_;}
    float metpt(){return metpt_;}
    float metphi(){return metphi_;}
    float mvacov00(){return mvacov00_;}
    float mvacov01(){return mvacov01_;}
    float mvacov10(){return mvacov10_;}
    float mvacov11(){return mvacov11_;}
};
typedef std::vector<Wmet> WmetCollection;


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

class Wjet{

    float pt_ = -999;
    float eta_ = -999;
    float phi_ = -999;
    float id_ = -999;
    float bptag_ = -999;
    float csvtag_ = -999;
    float bjet_ = -999;
    float jecfactor_ = -999;
    float jetlooseID_ = -999;
    float pujetetaid_ = -999;
  public:

    void pt(float x){pt_ = x;}
    void eta(float x){eta_ = x;}
    void phi(float x){phi_ = x;}
    void id(float x){id_ = x;}
    void bptag(float x){bptag_ = x;}
    void csvtag(float x){csvtag_ = x;}
    void bjet(float x){bjet_ = x;}
    void jecfactor(float x){jecfactor_ = x;}
    void jetlooseID(float x){jetlooseID_ = x;}
    void pujetetaid(float x){pujetetaid_ = x;}

    float pt(){return pt_;}
    float eta(){return eta_;}
    float phi(){return phi_;}
    float id(){return id_;}
    float bptag(){return bptag_;}
    float csvtag(){return csvtag_;}
    float bjet(){return bjet_;}
    float jecfactor(){return jecfactor_;}
    float jetlooseID(){return jetlooseID_;}
    float pujetetaid(){return pujetetaid_;}
};
typedef std::vector<Wjet> WjetCollection;


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

