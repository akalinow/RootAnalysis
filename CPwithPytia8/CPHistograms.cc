#include <iostream>
#include <cmath>

#include "CPHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"

CPHistograms::CPHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::CPHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::CPHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::~CPHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool CPHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill2DHistogram(name,val1,val2,weight)){
		if(name.find("Pt")!=std::string::npos) hTemplateName = "h2DPt";
		if(name.find("EtaHit")!=std::string::npos) hTemplateName = "h2DEtaHit";
		if(name.find("PhiHit")!=std::string::npos) hTemplateName = "h2DPhiHit";
		if(name.find("EtaVx")!=std::string::npos) hTemplateName = "h2DEtaVx";
		if(name.find("PhiVx")!=std::string::npos) hTemplateName = "h2DPhiVx";
		if(name.find("RateTot")!=std::string::npos) hTemplateName = "h2DRateTot";
		if(name.find("RateVsEta")!=std::string::npos) hTemplateName = "h2DRateVsEta";
		if(name.find("RateVsPt")!=std::string::npos) hTemplateName = "h2DRateVsPt";
		if(name.find("RateVsQuality")!=std::string::npos) hTemplateName = "h2DRateVsQuality";                                              
		std::cout<<"Adding histogram: "<<name<<std::endl;
		this->add2DHistogram(name,"",
				this->get2DHistogram(hTemplateName)->GetNbinsX(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
				this->get2DHistogram(hTemplateName)->GetNbinsY(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmax(),
				file_);
		return AnalysisHistograms::fill2DHistogram(name,val1,val2,weight);
	}
	return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
 //Make template histos
 add2DHistogram("h2DPt","",150,0,150,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaHit","",8*26,0.0,1.6,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHit","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaVx","",8*25,0.0,1.6,2,-0.5,1.6,file_);
 add2DHistogram("h2DPhiVx","",4*32,-0.2,2.2,2,-0.5,1.5,file_);
 //Rate histos
 add2DHistogram("h2DRateTot","",400,1,201,142,0,142,file_);
 add2DHistogram("h2DRateVsEta","",400,1,201,25,0.0,1.6,file_);

 add2DHistogram("h2DRateVsPt","",400,1,201,60,0,30,file_);

 add2DHistogram("h2DRateVsQuality","",400,1,201,10,-0.5,9.5,file_);

   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::finalizeHistograms(int nRuns, float weight){
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
