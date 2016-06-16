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
HTTWeightHistograms::HTTWeightHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::HTTWeightHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTWeightHistograms::~HTTWeightHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HTTWeightHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";
  if(name.find("h1DNPV")!=std::string::npos) templateName = "h1DNPVTemplate";
  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTWeightHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
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
