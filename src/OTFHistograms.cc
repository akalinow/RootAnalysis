#include <iostream>

#include "OTFHistograms.h"


int nPtBins = 32;
float ptBinsTmp[33]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
  		 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
  		 160. };

OTFHistograms::OTFHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::OTFHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::OTFHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::~OTFHistograms(){ 

  std::cout<<__func__<<std::endl;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

    add2DHistogram("h2DPt","",100,0,100,2,-0.5,1.5,file_);


   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(int nRuns, float weight){


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
