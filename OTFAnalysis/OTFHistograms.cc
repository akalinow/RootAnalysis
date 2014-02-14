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
    //Gmt
    add2DHistogram("h2DGmtPt","",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaHit","",230,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiHit","",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaVx","",230,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiVx","",64,-3.2,3.2,2,-0.5,1.5,file_);
    //Otf
    add2DHistogram("h2DOtfPt","",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaHit","",230,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiHit","",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaVx","",230,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiVx","",64,-3.2,3.2,2,-0.5,1.5,file_);

   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(int nRuns, float weight){


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
