//
// Original Author:  Artur Kalinowski
//         Created:  Tue Oct 24 15:08:51 CEST 2006
// $Id: AnalysisHistograms.cc,v 1.8 2010/09/14 11:37:23 cbern Exp $
//
//
// system include files
#include <memory>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "omp.h"
#include "AnalysisHistograms.h"

//
// constructors and destructor
//
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

std::vector<AnalysisHistograms*> AnalysisHistograms::container;
AnalysisHistograms::~AnalysisHistograms(){
 
  using namespace std;

  finalizeHistograms();

}


void AnalysisHistograms::addProfile(const std::string& name, 
				    const std::string& title, 
				    int nBinsX, float xlow, float xhigh, 
				    const TFileDirectory* myDir) {
  
  using namespace std;
  
  TProfile *hTmp = myDir->make<TProfile>(name.c_str(),title.c_str(),
					 nBinsX,xlow,xhigh);  
  if(myProfiles_.find(name)==myProfiles_.end()) myProfiles_[name] = hTmp;
  else cout<<"ERROR Substituting existing profile!"<<endl;
  
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title,
			                             int nBinsX, float xlow, float xhigh,
					                     const TFileDirectory* myDir){

  using namespace std;

  TH1F *hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);//myDir->make<TH1F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);
  hTmp->Sumw2();

  if(my1Dhistograms_.find(name)==my1Dhistograms_.end()) my1Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;

}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title, int nBinsX, float* bins,
					 const TFileDirectory* myDir){

  using namespace std;
  TH1F *hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,bins);//myDir->make<TH1F>(name.c_str(),title.c_str(),nBinsX,bins);
  hTmp->Sumw2();
  if(my1Dhistograms_.find(name)==my1Dhistograms_.end()) my1Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}

//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					int nBinsX, float xlow, float xhigh,
					 int nBinsY, float ylow, float yhigh,
					 const TFileDirectory* myDir){
   
  using namespace std;

  TH2F *hTmp= new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh, nBinsY, ylow,yhigh); 
  hTmp->SetDirectory(0);  
  hTmp->Sumw2();

  if(my2Dhistograms_.find(name)==my2Dhistograms_.end()) my2Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;

}

TH2F* AnalysisHistograms::registeredCopy(TH2F * h2f, std::string name, const TFileDirectory* myDir){
  
 TH2F* h2f1 =myDir->make<TH2F>(name.c_str(),"",
			      h2f->GetNbinsX(),h2f->GetXaxis()->GetXmin(),
				h2f->GetXaxis()->GetXmax(),h2f->GetNbinsY(),
				h2f->GetYaxis()->GetXmin(),h2f->GetYaxis()->GetXmax()); 
 h2f1->Sumw2();
 return h2f1;
  
}
TH3F* AnalysisHistograms::registeredCopy(TH3F * h2f, std::string name, const TFileDirectory* myDir){
  //Not sure about this, no data to test
 TH3F* h2f1 =myDir->make<TH3F>(name.c_str(),"",h2f->GetNbinsX(),h2f->GetXaxis()->GetXmin(),h2f->GetXaxis()->GetXmax(),
			       h2f->GetNbinsY(),h2f->GetYaxis()->GetXmin(),h2f->GetYaxis()->GetXmax(),
			       h2f->GetNbinsZ(),h2f->GetZaxis()->GetXmin(),h2f->GetZaxis()->GetXmax());

 h2f1->Sumw2();
 return h2f1;
  
}
TH1F* AnalysisHistograms::registeredCopy(TH1F * h2f, std::string name, const TFileDirectory* myDir){
  
 TH1F* h2f1 =myDir->make<TH1F>(name.c_str(),"",h2f->GetNbinsX(),h2f->GetXaxis()->GetXmin(),h2f->GetXaxis()->GetXmax());
 h2f1->Sumw2();
 return h2f1;
  
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float* binsX,
					 int nBinsY, float* binsY, 					 
					 const TFileDirectory* myDir){
  using namespace std;
  
  TH2F *hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);
  //myDir->make<TH2F>(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);
   hTmp->SetDirectory(0);  
  hTmp->Sumw2();

  if(my2Dhistograms_.find(name)==my2Dhistograms_.end()) my2Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float xlow, float xhigh,
					 int nBinsY, double* binsY, 					 
					 const TFileDirectory* myDir){

  using namespace std;
  
  TH2F *hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);
  //myDir->make<TH2F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);
   hTmp->SetDirectory(0);  
  hTmp->Sumw2();

   if(my2Dhistograms_.find(name)==my2Dhistograms_.end()) my2Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float xlow, float xhigh,
					 int nBinsY, float ylow, float yhigh,
					 int nBinsZ, float zlow, float zhigh,
					 const TFileDirectory* myDir){
  using namespace std;
  
  TH3F *hTmp =  new TH3F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);
  //myDir->make<TH3F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);
  hTmp->Sumw2();

  if(my3Dhistograms_.find(name)==my3Dhistograms_.end()) my3Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, double* binsX,
					 int nBinsY, double* binsY,
					 int nBinsZ, double* binsZ,
					 const TFileDirectory* myDir){

  using namespace std;
  
  TH3F *hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);
  //myDir->make<TH3F>(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);
  hTmp->Sumw2();

  if(my3Dhistograms_.find(name)==my3Dhistograms_.end()) my3Dhistograms_[name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}



void AnalysisHistograms::fillProfile(const std::string& name, float x, float val, float weight) {
  using namespace std;

  if(myProfiles_.find(name)!=myProfiles_.end()) 
    myProfiles_[name]->Fill(x, val,weight);
  else cout<<"ERROR: profile : "<<name<<" not found!"<<endl;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  AnalysisHistograms::fill1DHistogram(const std::string& name, float val, float weight){

 using namespace std;

  if(my1Dhistograms_.find(name)!=my1Dhistograms_.end()) my1Dhistograms_[name]->Fill(val,weight);
  else{
	 // cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
	  return false;
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

 using namespace std;

  if(my2Dhistograms_.find(name)!=my2Dhistograms_.end()) my2Dhistograms_[name]->Fill(val1,val2,weight);
  else{
	 // cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
	  return false;
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
bool  AnalysisHistograms::fill3DHistogram(const std::string& name, float val1, float val2, float val3,
					  float weight){

 using namespace std;

 if(my3Dhistograms_.find(name)!=my3Dhistograms_.end()) my3Dhistograms_[name]->Fill(val1,val2,val3,weight);
  else{
	 //cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
	  return false;
  }
 return true;
}
//////////////////////////////////////////////////////////////////////////////


TProfile* AnalysisHistograms::getProfile(const std::string& name) {

 using namespace std;

 if(myProfiles_.find(name)!=myProfiles_.end()) return (TProfile*)(myProfiles_[name]->Clone());
 else cout<<"ERROR: Profile : "<<name<<" not found!"<<endl;
 return 0;

}

TH1F* AnalysisHistograms::get1DHistogram(const std::string& name){

 using namespace std;

 if(my1Dhistograms_.find(name)!=my1Dhistograms_.end()) return (TH1F*)(my1Dhistograms_[name]->Clone());
// else cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
 return 0;

}
//////////////////////////////////////////////////////////////////////////////
TH2F* AnalysisHistograms::get2DHistogram(const std::string& name, bool noClone){

 using namespace std;

 if(noClone && my2Dhistograms_.find(name)!=my2Dhistograms_.end()) return my2Dhistograms_[name];
 else if(my2Dhistograms_.find(name)!=my2Dhistograms_.end()) return (TH2F*)(my2Dhistograms_[name]->Clone());
 //else cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
 return 0;

}



//////////////////////////////////////////////////////////////////////////////
TH3F* AnalysisHistograms::get3DHistogram(const std::string& name){

 using namespace std;

 if(my3Dhistograms_.find(name)!=my3Dhistograms_.end()) return (TH3F*)(my3Dhistograms_[name]->Clone());
 else cout<<"ERROR: Histogram: "<<name<<" not found!"<<endl;
 return 0;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::resetHistos(std::pair<const std::string, TH1*> aPair){
  aPair.second->Reset();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::init(TFileDirectory *myDir,
			      const std::string & name){ 

  name_ = name;
  myDirCopy = myDir;
  if(!histosInitialized_){
    if(name_.size()){
      mySecondaryDirs_.push_back(myDir->mkdir(name_));
      file_ = &mySecondaryDirs_[mySecondaryDirs_.size()-1];
    }
    else file_ = myDir;
    defineHistograms();
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::finalizeHistograms(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::clear(){

  for_each(myProfiles_.begin(), myProfiles_.end(), &AnalysisHistograms::resetHistos);
  for_each(my1Dhistograms_.begin(), my1Dhistograms_.end(), &AnalysisHistograms::resetHistos);
  for_each(my2Dhistograms_.begin(), my2Dhistograms_.end(), &AnalysisHistograms::resetHistos);
  for_each(my3Dhistograms_.begin(), my3Dhistograms_.end(), &AnalysisHistograms::resetHistos);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::write(){

  finalizeHistograms();
  //file->Write();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double* AnalysisHistograms::equalRanges(int nSteps, double min, double max, double *ranges){

  for(int i=0;i<=nSteps;++i){
    double val = min + i*(max-min)/nSteps;
    ranges[i] = val;
  }

  return ranges;
}

template <typename THmap> void AnalysisHistograms:: pieceHistogramsTogether(THmap &map, std::vector<THmap*> & array){
  for(int i =0; i < array.size();++i){
    for(auto element: *array[i]){
     if(map.find(element.first) == map.end())
	map[element.first]= registeredCopy(element.second, element.first,myDirCopy);
     for(int a =0; a< map[element.first]->GetSize();++a){
		(*map[element.first])[a] += (*element.second)[a];
		 map[element.first]->SetBinError(a,  map[element.first]->GetBinError(a) +
		  element.second->GetBinError(a) * element.second->GetBinError(a) );
      }
      delete element.second;
    } 
    array[i]->clear();
   } 
}
template <typename THmap> void AnalysisHistograms::clearTemplatesFromMaps( THmap & map, std::vector<THmap*> & array){
  for(int i =0; i < array.size();++i)
     for(auto &element:map)
       array[i]->erase(element.first);
}
template <typename THmap> void AnalysisHistograms::fixBinErrors( THmap & map){
   for(auto &element: map){
      for(int a =0; a< element.second->GetSize();++a){
	element.second->SetBinError(a, sqrt(element.second->GetBinError(a)));
      }
   }
}

  
 void AnalysisHistograms::pieceHistogramsTogether2D(std::vector<std::map<std::string,TH2F*>*> & map)
{
  clearTemplatesFromMaps<std::map<std::string,TH2F*>>(my2Dhistograms_,map);
  pieceHistogramsTogether<std::map<std::string,TH2F*>>(my2Dhistograms_, map);
  fixBinErrors<std::map<std::string,TH2F*>>(my2Dhistograms_);
}
 void AnalysisHistograms::pieceHistogramsTogether1D(std::vector<std::map<std::string,TH1F*>*> & map)
{
  clearTemplatesFromMaps<std::map<std::string,TH1F*>>(my1Dhistograms_,map);
  pieceHistogramsTogether<std::map<std::string,TH1F*>>(my1Dhistograms_, map);
  fixBinErrors<std::map<std::string,TH1F*>>(my1Dhistograms_);
}
 void AnalysisHistograms::pieceHistogramsTogether3D(std::vector<std::map<std::string,TH3F*>*> & map)
{
  clearTemplatesFromMaps<std::map<std::string,TH3F*>>(my3Dhistograms_,map);
  pieceHistogramsTogether<std::map<std::string,TH3F*>>(my3Dhistograms_, map);
  fixBinErrors<std::map<std::string,TH3F*>>(my3Dhistograms_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
