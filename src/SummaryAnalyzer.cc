#include <sstream>
#include <algorithm>
#include <string>
#include <bitset>
#include <vector>
#include <cmath>

#include "SummaryAnalyzer.h"
#include "TreeAnalyzer.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TObjArray.h"
#include "TCollection.h"
#include "TBits.h"
#include "TObject.h"
#include "TROOT.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
SummaryAnalyzer::SummaryAnalyzer(const std::string & aName):Analyzer(aName){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
SummaryAnalyzer::~SummaryAnalyzer(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* SummaryAnalyzer::clone() const{

  SummaryAnalyzer* clone = new SummaryAnalyzer(name());
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SummaryAnalyzer::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  myDir_ = &aDir;

  mySelections_ = aSelections;

  eventWeight_ = 1.0;

  mySelectionsTree_ = aDir.make<TTree>("tree","Selections bit words");
  branchWeight_ = mySelectionsTree_->Branch("eventWeight", &eventWeight_);

  //selectionWord_ = TBits(mySelections_->bits().size());
  //bitsBranch_ = mySelectionsTree_->Branch("selectionWord",&selectionWord_);

  /////////////////////////////////////////
  ///Create histogram with the cuts names
  int nCuts = (int)(mySelections_->bits().size());
  int nBins = nCuts+1;

  hCutNames_ = aDir.make<TH1F>("hCutNames","Names of the selections",nBins,-0.5,nBins-0.5);

  const std::vector<std::string> cutNames = mySelections_->strings();
  std::vector<std::string>::const_iterator ci = cutNames.begin();
  int iBinX = 1;
  hCutNames_->GetXaxis()->SetBinLabel(iBinX,"preselection");
  for(;ci!=cutNames.end();++ci){
    ++iBinX;
    hCutNames_->GetXaxis()->SetBinLabel(iBinX,ci->c_str());
  }
  ///////////////////////////////////////////////////////
  /////////Create the cut counters for all the cut flow flavours
  for(unsigned int i=0;i<selectionFlavours_.size();++i){
    cutCounters_.push_back(float(1.0));	  
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SummaryAnalyzer::finalize(){ 

  using namespace std;

  for(unsigned int i=0;i<selectionFlavours_.size();++i){
      fillEffHisto(selectionFlavours_[i]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  SummaryAnalyzer::fillEffHisto(std::string type){

  std::string label = "";

  int nCuts = (int)(mySelections_->bits().size());

  const std::vector<std::string> & cutNames = mySelections_->strings();
  std::vector<std::string>::const_iterator ci = cutNames.begin();
  //////////////////
  std::vector<TBits> typeBits;
  TBits tmpBits(nCuts);
  for(int i=0;i<nCuts;++i){    
    if(ci->find(type)==std::string::npos){
      ++ci;
      continue;
    }    
    tmpBits.SetBitNumber(i);
    ++ci;
    typeBits.push_back(*(TBits*)tmpBits.Clone());
  }
  TObjArray* branchList = mySelectionsTree_->GetListOfBranches();
  int maxIndex = branchList->GetLast();
  for(int i=2;i<=maxIndex;++i){    
    std::string branchName(branchList->At(i)->GetName());
    //cout<<"branch "<<branchName<<endl;
    if(branchName.find(type)==std::string::npos) continue;
    float tmpVal = 0;
    TBranch *tmpBranch = (TBranch*)branchList->At(i);
    tmpBranch->SetAddress(&tmpVal);
    std::string hName = "h"+branchName;
    ////
    TH1F *h1D = (TH1F*)cutCounterHistos_.FindObject(hName.c_str());
    TH2F *h = 0; 
    if(h1D && h1D->GetDimension()==1){
    h1D->SetDirectory(0);
    int nBins = nCuts+1;
    h = myDir_->make<TH2F>(hName.c_str(),
			   h1D->GetTitle(),h1D->GetNbinsX(),
			   h1D->GetXaxis()->GetXmin(),
			   h1D->GetXaxis()->GetXmax(),  
			   nBins, -0.5,nBins-0.5);
    h->SetXTitle(h1D->GetXaxis()->GetTitle());   
    h->Sumw2();
    }
    else{
      //std::cout<<"ERROR Histogram: "<<hName<<" not found!"<<std::endl;
      continue;
    }
    for(int k=0;k<tmpBranch->GetEntries();++k){
      tmpBranch->GetEntry(k);
      bitsBranch_->GetEntry(k);
      branchWeight_->GetEntry(k);
      int valY=0;
      for(uint l=0;l<typeBits.size();++l) valY+= ((selectionWord_ & typeBits[l])==typeBits[l]);      
      h->Fill(tmpVal,valY,eventWeight_);
    }

    /////////////////////
    for(int iBinX=0;iBinX<=h->GetNbinsX()+1;iBinX++){
      TH1D *hTmp = h->ProjectionY("hTmp",iBinX,iBinX);
      hTmp->SetDirectory(0);
      TH1D *hIntegrated = Integrate(hTmp);
      hIntegrated->SetDirectory(0);
      for(int iBinY=0;iBinY<=hIntegrated->GetNbinsX()+1;iBinY++){
	h->SetBinContent(iBinX,iBinY,hIntegrated->GetBinContent(iBinY));
	h->SetBinError(iBinX,iBinY,hIntegrated->GetBinError(iBinY));
      }
      delete hTmp;
      delete hIntegrated;
    }
    ci = cutNames.begin();
    for(int i=0;i<nCuts;++i){    
      if(ci->find(type)==std::string::npos){
	++ci;
	continue;
      }    
      h->GetYaxis()->SetBinLabel(i+2,ci->c_str());
      ++ci;
    }
    h->GetYaxis()->SetBinLabel(1,"preselection");
    clean2DHisto(h);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool SummaryAnalyzer::analyze(const EventProxyBase& iEvent){
  /*
  const pat::strbitset::bit_vector& orderedBits = mySelections_->bits(); 
  pat::strbitset::bit_vector::const_iterator ci = orderedBits.begin();
  
  int counter = 0;
  selectionWord_.ResetAllBits();
  for(;ci!=orderedBits.end();++ci) selectionWord_.SetBitNumber(counter++,*ci);
  */  
  ////
  mySelectionsTree_->Fill();

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SummaryAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
   tree->Branch(("CutCounter"+selectionFlavours_[i]).c_str(),&cutCounters_[i]);
 }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SummaryAnalyzer::addCutHistos(TList *aList){

  for(unsigned int i=0;i<selectionFlavours_.size();++i){
    std::string hName = "hCutCounter"+selectionFlavours_[i];
    std::string hTitle = "Cut counter for "+selectionFlavours_[i]+" selections; ; #sigma [fb]";
    cutCounterHistos_.Add(new TH1F(hName.c_str(),hTitle.c_str(),1,0.5,1.5));
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Clean non labeled bins
void SummaryAnalyzer::cleanHisto(TH2F *hCutCounter){

  int nBinsFill = 0;
  for(int i=0;i<hCutCounter->GetNbinsY()+1;++i){
    TString tmpStr(hCutCounter->GetYaxis()->GetBinLabel(i));
    if(tmpStr.Length()>1) nBinsFill++;
  }
  std::string tmpStr(hCutCounter->GetName());
  tmpStr+="1D";
  TH1F *hCutCounterTmp = myDir_->make<TH1F>(tmpStr.c_str(),
					    hCutCounter->GetTitle(),
					    nBinsFill, -0.5,nBinsFill-0.5);
  int iBin = 0;
  for(int i=0;i<hCutCounter->GetNbinsY()+1;++i){
    TString tmpStr(hCutCounter->GetYaxis()->GetBinLabel(i));
    if(tmpStr.Length()>1){
       iBin++;
      hCutCounterTmp->GetXaxis()->SetBinLabel(iBin,
					      hCutCounter->GetYaxis()->GetBinLabel(i));
     }
     if(i<=nBinsFill){
	 hCutCounterTmp->SetBinContent(i,hCutCounter->GetBinContent(1,i));
	 hCutCounterTmp->SetBinError(i,hCutCounter->GetBinError(1,i));
     }
  }
  //////////////////////////
  hCutCounter->SetDirectory(0);
  delete hCutCounter;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Clean non labeled bins
void SummaryAnalyzer::clean2DHisto(TH2F *hCutCounter){

  std::string hName(hCutCounter->GetName());

  int nBinsFill = 0;
  for(int i=0;i<hCutCounter->GetNbinsY()+1;++i){
    TString tmpStr(hCutCounter->GetYaxis()->GetBinLabel(i));
    if(tmpStr.Length()>1) nBinsFill++;
  }
  int nBinsX = hCutCounter->GetNbinsX();
  float lowX = hCutCounter->GetXaxis()->GetXmin();
  float highX = hCutCounter->GetXaxis()->GetXmax();

  hCutCounter->SetDirectory(0);
  TH2F *hCutCounterTmp = myDir_->make<TH2F>((hName+"Clean").c_str(),hCutCounter->GetTitle(),
					    nBinsX, lowX, highX,
					    nBinsFill, -0.5,nBinsFill-0.5);  
  int iBinY = 0;
  for(int i=0;i<hCutCounter->GetNbinsY()+1;++i){
    TString tmpStr(hCutCounter->GetYaxis()->GetBinLabel(i));
     if(tmpStr.Length()>1){
      iBinY++;     
	hCutCounterTmp->GetYaxis()->SetBinLabel(iBinY,
						hCutCounter->GetYaxis()->GetBinLabel(i));
     }
     if(i<=nBinsFill){
       for(int iBinX=0;iBinX<=hCutCounter->GetNbinsX()+1;++iBinX){
	 hCutCounterTmp->SetBinContent(iBinX,i,hCutCounter->GetBinContent(iBinX,i));
	 hCutCounterTmp->SetBinError(iBinX,i,hCutCounter->GetBinError(iBinX,i));
       }
     }
  }
  //hCutCounter->Print("all");
  //hCutCounterTmp->Print("all");
  //////////////////////////
  if(hName.find("CutCounter")!=std::string::npos){
    cleanHisto(hCutCounter);
    hCutCounterTmp->SetDirectory(0);
    delete hCutCounterTmp;
  }
  else delete hCutCounter;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TH1D * SummaryAnalyzer::Integrate(TH1D * histoD) {

   TH1D * histoI = new TH1D(*histoD); 
   Float_t *  cont = new Float_t [histoD->GetNbinsX()+2];  //with under+overflow
   Float_t *  errs = new Float_t [histoD->GetNbinsX()+2];  //with under+overflow
   histoI->Reset();
   
// bin=0 underf
// bin 1-GetNbinsX() -conten
// bin GetNbinsX()+1 overflow

   Int_t i;
   for (i = 0; i <= histoD->GetNbinsX()+1; i++) { 
      cont[i] = histoD->GetBinContent(i);   
      errs[i] = histoD->GetBinError(i);
   }
   Float_t sum=0.;
   Float_t sume2=0.;
   ////////////////////
   for (i = histoD->GetNbinsX()+1; i >= 0; i--) {
        sum+=cont[i];
        sume2+=errs[i]*errs[i];
	///////////////////////////
        histoI->SetBinContent(i,sum);
        histoI->SetBinError(i,sqrt(sume2));
   }   
   return histoI;
}
////////////////////////////////////////////////////
