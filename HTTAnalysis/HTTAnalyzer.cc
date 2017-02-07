#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "MuTauSpecifics.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName, const std::string & aDecayMode) : Analyzer(aName){

        //pileupCalc.py -i lumiSummary_Run2016BCDE_PromptReco_v12.json
        //--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
        //--calcMode true --minBiasXsec 69200 --maxPileupBin 60 --numPileupBins 600 Data_Pileup_Cert_271036-277148.root

        ///Load ROOT file with PU histograms.
        std::string filePath = "Data_Pileup_2016_BCDEFG_v26.root";
        //filePath = "Data_Pileup_2016_July22.root";
        puDataFile_ = new TFile(filePath.c_str());

        filePath = "MC_Spring16_PU25ns_V1.root";
        puMCFile_ = new TFile(filePath.c_str());

        categoryDecisions.resize((int)HTTAnalyzer::DUMMY);

        myChannelSpecifics = new MuTauSpecifics(this);

        nPCAMin_ = 0.003;

        ntupleFile_ = 0;
        hStatsFromFile = 0;

        h2DMuonIdCorrections = 0;
        h2DMuonIsoCorrections = 0;
        h2DMuonTrgCorrections = 0;
        h3DTauCorrections = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

        if(myHistos_) delete myHistos_;
        if(puDataFile_) delete puDataFile_;
        if(puMCFile_) delete puMCFile_;

        if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
        if(h2DMuonIsoCorrections) delete h2DMuonIsoCorrections;
        if(h2DMuonTrgCorrections) delete h2DMuonTrgCorrections;
        if(h3DTauCorrections) delete h3DTauCorrections;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTAnalyzer::clone() const {

        std::string myDecayMode = myChannelSpecifics->getDecayModeName();
        HTTAnalyzer* clone = new HTTAnalyzer(name(),myDecayMode);
        clone->setHistos(myHistos_);
        return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initialize(TDirectory* aDir,
                             pat::strbitset *aSelections){

        mySelections_ = aSelections;

        myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initializeCorrections(){

#pragma omp critical
        {
                std::string filePath = "htt_scalefactors_v5.root";
                TFile aFile(filePath.c_str());
                RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile.Get("w");

                RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
                RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio");
                RooAbsReal *muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data");//OR of the HLT_IsoMu22 and HLT_IsoTkMu22
                //RooAbsReal *tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_m_pt30_sf");

                h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
                                                                                   *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                   RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                   RooFit::Scaling(kFALSE));

                h2DMuonIsoCorrections = (TH2F*)muon_iso_scalefactor->createHistogram("h2DMuonIsoCorrections",
                                                                                     *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                     RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                     RooFit::Scaling(kFALSE));

                h2DMuonTrgCorrections = (TH2F*)muon_trg_efficiency->createHistogram("h2DMuonTrgCorrections",
                                                                                    *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                    RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                    RooFit::Scaling(kFALSE));
        }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){
        myHistos_->finalizeHistograms(0,1.0);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> HTTAnalyzer::getSeparatedJets(const EventProxyHTT & myEventProxy,
                                                       float deltaR){

        std::vector<HTTParticle> separatedJets;

        for(auto aJet: *myEventProxy.jets) {
                float dRLeg2 = aJet.getP4().DeltaR(aLeg2.getP4());
                float dRLeg1 = aJet.getP4().DeltaR(aLeg1.getP4());
                bool loosePFJetID = aJet.getProperty(PropertyEnum::PFjetID)>=1;
                if(dRLeg1>deltaR && dRLeg2>deltaR && loosePFJetID) separatedJets.push_back(aJet);
        }

        return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

        aEvent = *myEventProxy.event;
        aPair = (*myEventProxy.pairs)[0];

        TLorentzVector met4v(aPair.getMET().X(),
                             aPair.getMET().Y(),
                             0,
                             aPair.getMET().Mod());

        aMET = HTTParticle();
        aMET.setP4(met4v);

        myChannelSpecifics->setAnalysisObjects(myEventProxy);

        aSeparatedJets = getSeparatedJets(myEventProxy, 0.5);
        aJet1 = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
        aJet2 = aSeparatedJets.size()>1 ? aSeparatedJets[1] : HTTParticle();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::addBranch(TTree *tree){/*tree->Branch("muonPt",&muonPt);*/
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillControlHistos(const std::string & hNameSuffix, float eventWeight,
                                    const sysEffects::sysEffectsEnum & aSystEffect){

        ///Histograms filled for each systematic effect
        ///Fill SVfit and visible masses
        const TLorentzVector & aVisSum = aLeg1.getP4(aSystEffect) + aLeg2.getP4(aSystEffect);
        float visMass = aVisSum.M();
        float higgsPt =  (aVisSum + aMET.getP4(aSystEffect)).Pt();
        float jetsMass = 0;
        if(nJets30>1) jetsMass = (aJet1.getP4(aSystEffect)+aJet2.getP4(aSystEffect)).M();

        myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4(aSystEffect).M(),eventWeight);
        myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix, visMass, eventWeight);
        myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(aSystEffect),eventWeight);

        ///Unrolled distributions for 2D fit
        myHistos_->fill2DUnrolledHistogram("h1DUnRollTauPtMassVis"+hNameSuffix, visMass, aLeg2.getP4(aSystEffect).Pt(),eventWeight);
        myHistos_->fill2DUnrolledHistogram("h1DUnRollHiggsPtMassSV"+hNameSuffix, aPair.getP4(aSystEffect).M(), higgsPt, eventWeight);
        myHistos_->fill2DUnrolledHistogram("h1DUnRollMjjMassSV"+hNameSuffix, aPair.getP4(aSystEffect).M(), jetsMass, eventWeight);

        myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aLeg1.getProperty(PropertyEnum::combreliso),eventWeight);
        if(aSystEffect!=sysEffects::NOMINAL_SVFIT) return;

        fillDecayPlaneAngle(hNameSuffix, eventWeight, aSystEffect);
        fillGenDecayPlaneAngle(hNameSuffix+"_Gen", eventWeight);
        fillDecayPlaneAngle(hNameSuffix+"_RefitPV", eventWeight);
        fillDecayPlaneAngle(hNameSuffix+"_AODPV", eventWeight);
        fillDecayPlaneAngle(hNameSuffix+"_GenPV", eventWeight);
        fillVertices(hNameSuffix+"_RefitPV", eventWeight);
        fillVertices(hNameSuffix+"_AODPV", eventWeight);

        ///Fill histograms with number of PV.
        myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.getNPV(),eventWeight);

        ///Fill leg1 (muon or first tau)
        myHistos_->fill1DHistogram("h1DPtLeg1"+hNameSuffix,aLeg1.getP4(aSystEffect).Pt(),eventWeight);
        myHistos_->fill1DHistogram("h1DEtaLeg1"+hNameSuffix,aLeg1.getP4(aSystEffect).Eta(),eventWeight);
        myHistos_->fill1DHistogram("h1DPhiLeg1"+hNameSuffix,aLeg1.getP4(aSystEffect).Phi(),eventWeight);
        myHistos_->fill1DHistogram("h1DIsoLeg1"+hNameSuffix,aLeg1.getProperty(PropertyEnum::combreliso),eventWeight);
        myHistos_->fill1DHistogram("h1DIDLeg1"+hNameSuffix,aLeg1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
        myHistos_->fill1DHistogram("h1DPtLeg1LeadingTk"+hNameSuffix,aLeg1.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);
        myHistos_->fill1DHistogram("h1DStatsLeg1DecayMode"+hNameSuffix, aLeg1.getProperty(PropertyEnum::decayMode), eventWeight);
        myHistos_->fill1DHistogram("h1DnPCALeg1"+hNameSuffix,aLeg1.getPCARefitPV().Mag(),eventWeight);

        ///Fill leg2 (tau)
        myHistos_->fill1DHistogram("h1DPtLeg2"+hNameSuffix,aLeg2.getP4(aSystEffect).Pt(),eventWeight);
        myHistos_->fill1DHistogram("h1DEtaLeg2"+hNameSuffix,aLeg2.getP4(aSystEffect).Eta(),eventWeight);
        myHistos_->fill1DHistogram("h1DPhiLeg2"+hNameSuffix,aLeg2.getP4(aSystEffect).Phi(),eventWeight);
        myHistos_->fill1DHistogram("h1DIDLeg2"+hNameSuffix,aLeg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
        myHistos_->fill1DHistogram("h1DStatsLeg2DecayMode"+hNameSuffix, aLeg2.getProperty(PropertyEnum::decayMode), eventWeight);
        myHistos_->fill1DHistogram("h1DnPCALeg2"+hNameSuffix,aLeg2.getPCARefitPV().Mag(),eventWeight);
        myHistos_->fill1DHistogram("h1DPtLeg2LeadingTk"+hNameSuffix,aLeg2.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);
        myHistos_->fill1DHistogram("h1DPtLeg1Leg2MET"+hNameSuffix,higgsPt,eventWeight);

        ///Fill jets info
        myHistos_->fill1DHistogram("h1DStatsNJ30"+hNameSuffix,nJets30,eventWeight);
        if(nJets30>0) {
                myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet1.getP4(aSystEffect).Pt(),eventWeight);
                myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet1.getP4(aSystEffect).Eta(),eventWeight);
                myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet1.getProperty(PropertyEnum::bCSVscore),eventWeight);
        }
        if(nJets30>1) myHistos_->fill1DHistogram("h1DWideMass2J"+hNameSuffix,jetsMass,eventWeight);
        myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aMET.getP4(aSystEffect).Pt(),eventWeight);

        if(aJet1.getProperty(PropertyEnum::bCSVscore)>0.8) {
                myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet1.getP4(aSystEffect).Pt(),eventWeight);
                myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet1.getP4(aSystEffect).Eta(),eventWeight);
        }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::fillVertices(const std::string & sysType, float eventWeight){

        TVector3 aVertexGen = aEvent.getGenPV();
        TVector3 aVertex = aEvent.getRefittedPV();
        if(sysType.find("AODPV")!=std::string::npos) aVertex = aEvent.getAODPV();
        if(sysType.find("RefitPV")!=std::string::npos) aVertex = aEvent.getRefittedPV();

        float pullX = aVertexGen.X() - aVertex.X();
        float pullY = aVertexGen.Y() - aVertex.Y();
        float pullZ = aVertexGen.Z() - aVertex.Z();

        myHistos_->fill1DHistogram("h1DVxPullX"+sysType,pullX, eventWeight);
        myHistos_->fill1DHistogram("h1DVxPullY"+sysType,pullY, eventWeight);
        myHistos_->fill1DHistogram("h1DVxPullZ"+sysType,pullZ, eventWeight);

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::passCategory(const HTTAnalyzer::muTauCategory & aCategory){

        if(categoryDecisions.size()==0) return false;
        else return categoryDecisions[(int)aCategory];

        return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

        const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
        sampleName = getSampleName(myEventProxy);

        std::string hNameSuffix = sampleName;
        float puWeight = getPUWeight(myEventProxy);
        float genWeight = getGenWeight(myEventProxy);
        float ptReweight = 1.0;
        if(sampleName.find("DYJets")!=std::string::npos ||
           sampleName.find("TTbar")!=std::string::npos)
                ptReweight = myEventProxy.event->getPtReWeight();
        float eventWeight = puWeight*genWeight*ptReweight;

        //Fill bookkeeping histogram. Bin 1 holds sum of weights.
        myHistos_->fill1DHistogram("h1DStats"+sampleName,0,eventWeight);
        myHistos_->fill1DHistogram("h1DNPartons"+hNameSuffix,myEventProxy.event->getLHEnOutPartons(),eventWeight);
        getPreselectionEff(myEventProxy);

        if(!myEventProxy.pairs->size()) return true;

        setAnalysisObjects(myEventProxy);

        std::pair<bool, bool> goodDecayModes = myChannelSpecifics->checkTauDecayMode(myEventProxy);
        bool goodGenDecayMode = goodDecayModes.first;
        bool goodRecoDecayMode = goodDecayModes.second;

        if(goodGenDecayMode) fillGenDecayPlaneAngle(sampleName+"GenNoOfflineSel", eventWeight);

        int tauIDmask = 0;
        for(unsigned int iBit=0; iBit<aEvent.ntauIds; iBit++) {
                if(aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
                if(aEvent.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
                if(aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }
        bool tauID = ( (int)aLeg2.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool muonKinematics = aLeg1.getP4().Pt()>24 && fabs(aLeg1.getP4().Eta())<2.1;

        bool trigger = aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
                       aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
                       aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
                       aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1);

        if(sampleName!="Data") trigger = true; //MC trigger included in muon SF
        if(!muonKinematics || !tauID || !trigger) return true;

        bool SS = aLeg2.getCharge()*aLeg1.getCharge() == 1;
        bool OS = aLeg2.getCharge()*aLeg1.getCharge() == -1;

        std::string categorySuffix = "";
        std::string systEffectName = "";
        for(unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT;
            iSystEffect<(unsigned int)sysEffects::DUMMY; ++iSystEffect) {

                sysEffects::sysEffectsEnum aSystEffect = static_cast<sysEffects::sysEffectsEnum>(iSystEffect);

                float leg1ScaleFactor = getLeptonCorrection(aLeg1.getP4(aSystEffect).Eta(), aLeg1.getP4(aSystEffect).Pt(), hadronicTauDecayModes::tauDecayMuon);
                float leg2ScaleFactor = getLeptonCorrection(aLeg2.getP4(aSystEffect).Eta(), aLeg2.getP4(aSystEffect).Pt(),
                                                           static_cast<hadronicTauDecayModes>(aLeg2.getProperty(PropertyEnum::decayMode)));
                float weightSyst = getSystWeight(aSystEffect);
                float eventWeightWithSyst=eventWeight*weightSyst*leg1ScaleFactor*leg2ScaleFactor;

                TLorentzVector met4v(aPair.getMET(aSystEffect).X(),
                                     aPair.getMET(aSystEffect).Y(), 0,
                                     aPair.getMET(aSystEffect).Mod());
                aMET.setP4(met4v);

                myChannelSpecifics->testAllCategories(aSystEffect);

                for(unsigned int iCategory = HTTAnalyzer::jet0; iCategory<HTTAnalyzer::W; ++iCategory) {
                        HTTAnalyzer::muTauCategory categoryType = static_cast<HTTAnalyzer::muTauCategory>(iCategory);

                        if(!passCategory(categoryType)) continue;

                        categorySuffix = std::to_string(iCategory);
                        systEffectName = HTTAnalyzer::systEffectName(iSystEffect);

                        if(systEffectName.find("CAT")!=std::string::npos) {
                                std::string categoryName = HTTAnalyzer::categoryName(iCategory);
                                systEffectName.replace(systEffectName.find("CAT"),3,categoryName);
                        }

                        if(OS) hNameSuffix = sampleName+"_OS_"+categorySuffix+systEffectName;
                        else if(SS) hNameSuffix = sampleName+"_SS_"+categorySuffix+systEffectName;

                        fillControlHistos(hNameSuffix, eventWeightWithSyst, aSystEffect);
                }
        }
        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
