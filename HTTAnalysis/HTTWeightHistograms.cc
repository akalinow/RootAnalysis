#include <iostream>
#include <cmath>

#include "HTTWeightHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"

HTTWeightHistograms::HTTWeightHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::~HTTWeightHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool HTTWeightHistograms::fill1DHistogram(const std::string& name, float val, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DNPV")!=std::string::npos) hTemplateName = "h1DNPVTemplate";
    std::cout<<"Adding histogram: "<<name<<" "<<file_<<" "<<file_->fullPath()<<std::endl;
    this->add1DHistogram(name,"",
			 this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmin(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmax(),
			 file_);
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   std::cout<<"Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::finalizeHistograms(int nRuns, float weight){

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
