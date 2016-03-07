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
AnalysisHistograms::~AnalysisHistograms(){
 
  using namespace std;

  //finalizeHistograms();

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void AnalysisHistograms::addProfile(const std::string& name, 
				    const std::string& title, 
				    int nBinsX, float xlow, float xhigh, 
				    const TFileDirectory* myDir) {
  
  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  
  TProfile *hTmp = 0;
  hTmp = new TProfile(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);
  //if(iThread==0) hTmp = myDir->make<TProfile>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);
  //else hTmp = new TProfile(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);

  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(myProfiles_[iThread].find(name)==myProfiles_[iThread].end()) myProfiles_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing profile!"<<endl;    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title,
			                             int nBinsX, float xlow, float xhigh,
					                     TFileDirectory* myDir){

  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  TH1F *hTmp = 0;
  hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);
  //if(iThread==0) hTmp = myDir->make<TH1F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);
  //else hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my1Dhistograms_[iThread].find(name)==my1Dhistograms_[iThread].end()) my1Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;  
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title, int nBinsX, float* bins,
					 const TFileDirectory* myDir){

  using namespace std;
  unsigned int iThread = omp_get_thread_num();

  TH1F *hTmp = 0;
  hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,bins);
  //if(iThread==0) hTmp = myDir->make<TH1F>(name.c_str(),title.c_str(),nBinsX,bins);
  //else hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,bins);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my1Dhistograms_[iThread].find(name)==my1Dhistograms_[iThread].end()) my1Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;  
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					int nBinsX, float xlow, float xhigh,
					 int nBinsY, float ylow, float yhigh,
					 const TFileDirectory* myDir){
   
  using namespace std;
  unsigned int iThread = omp_get_thread_num();

  TH2F *hTmp = 0;
  hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh, nBinsY, ylow,yhigh);
  //if(iThread==0) hTmp = myDir->make<TH2F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh, nBinsY, ylow,yhigh);
  //else hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh, nBinsY, ylow,yhigh);
  
  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float* binsX,
					 int nBinsY, float* binsY, 					 
					 const TFileDirectory* myDir){
  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  
  TH2F *hTmp =0;
  hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);
  //if(iThread==0) hTmp = myDir->make<TH2F>(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);
  //else hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float xlow, float xhigh,
					 int nBinsY, double* binsY, 					 
					 const TFileDirectory* myDir){

  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  
  TH2F *hTmp = 0;
  hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);
  //if(iThread==0) hTmp = myDir->make<TH2F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);
  //else hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, float xlow, float xhigh,
					 int nBinsY, float ylow, float yhigh,
					 int nBinsZ, float zlow, float zhigh,
					 const TFileDirectory* myDir){
  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  
  TH3F *hTmp = 0;
  hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);
  //if(iThread==0) hTmp = myDir->make<TH3F>(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);
  //else hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my3Dhistograms_[iThread].find(name)==my3Dhistograms_[iThread].end()) my3Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;  
}
//////////////////////////////////////////////////////////////////////////////
void  AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
					 int nBinsX, double* binsX,
					 int nBinsY, double* binsY,
					 int nBinsZ, double* binsZ,
					 const TFileDirectory* myDir){

  using namespace std;
  unsigned int iThread = omp_get_thread_num();
  
  TH3F *hTmp = 0;
  hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);
  //if(iThread==0) hTmp = myDir->make<TH3F>(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);
  //else hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);

  hTmp->Sumw2();
  //hTmp->SetDirectory(0);

  //if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir->getBareDirectory());
  if(my3Dhistograms_[iThread].find(name)==my3Dhistograms_[iThread].end()) my3Dhistograms_[iThread][name] = hTmp;
  else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fillProfile(const std::string& name, float x, float val, float weight) {
  using namespace std;

  unsigned int iThread = omp_get_thread_num();

  std::unordered_map<std::string,TProfile*>::iterator it = myProfiles_[iThread].find(name); 
  if(it!=myProfiles_[iThread].end()) it->second->Fill(x,val,weight);
  else return false;

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  AnalysisHistograms::fill1DHistogram(const std::string& name, float val, float weight){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();
 
 std::unordered_map<std::string,TH1F*>::iterator it = my1Dhistograms_[iThread].find(name);
 
 if(it!=my1Dhistograms_[iThread].end()) it->second->Fill(val,weight); 
 else return false;
 
 return true;
}
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();

 std::unordered_map<std::string,TH2F*>::iterator it = my2Dhistograms_[iThread].find(name); 
 if(it!=my2Dhistograms_[iThread].end()) it->second->Fill(val1,val2,weight);
 else return false;
  
 return true;
}
//////////////////////////////////////////////////////////////////////////////
bool  AnalysisHistograms::fill3DHistogram(const std::string& name, float val1, float val2, float val3,
					  float weight){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();

 std::unordered_map<std::string,TH3F*>::iterator it = my3Dhistograms_[iThread].find(name); 
 if(it!=my3Dhistograms_[iThread].end()) it->second->Fill(val1,val2,val3,weight);
 else return false;
   
 return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TProfile* AnalysisHistograms::getProfile(const std::string& name, bool noClone) {

 using namespace std;


 unsigned int iThread = omp_get_thread_num();
 if(name.find("Template")!=std::string::npos) iThread = 0;

 std::unordered_map<std::string,TProfile*>::iterator it = myProfiles_[iThread].find(name);

 if(noClone && it!=myProfiles_[iThread].end()) return it->second;
 else if(it!=myProfiles_[iThread].end()){
     TProfile* hClone =  (TProfile*)(it->second->Clone());
     hClone->SetDirectory(0);
     return hClone;
   }
 
 else cout<<"ERROR: Profile : "<<name<<" not found!"<<endl;
 return 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TH1F* AnalysisHistograms::get1DHistogram(const std::string& name, bool noClone){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();
 if(name.find("Template")!=std::string::npos) iThread = 0;

 std::unordered_map<std::string,TH1F*>::iterator it = my1Dhistograms_[iThread].find(name);
 
 if(noClone && it!=my1Dhistograms_[iThread].end()) return it->second;
 else if(it!=my1Dhistograms_[iThread].end()){
     TH1F* hClone =  (TH1F*)(it->second->Clone());
     hClone->SetDirectory(0);
     return hClone;
   }
 else cout<<"ERROR: Histogram: "<<name<<" not found in thread: "<<iThread<<endl;

 return 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TH2F* AnalysisHistograms::get2DHistogram(const std::string& name, bool noClone){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();
 if(name.find("Template")!=std::string::npos) iThread = 0;

 std::unordered_map<std::string,TH2F*>::iterator it = my2Dhistograms_[iThread].find(name);
 
 if(noClone && it!=my2Dhistograms_[iThread].end()) return it->second;
 else if(it!=my2Dhistograms_[iThread].end()){
   TH2F* hClone =  (TH2F*)(it->second->Clone());
   hClone->SetDirectory(0);
   return hClone;
 }
 else cout<<"ERROR: Histogram: "<<name<<" not found in thread: "<<iThread<<endl;
 return 0;
}
//////////////////////////////////////////////////////////////////////////////
TH3F* AnalysisHistograms::get3DHistogram(const std::string& name, bool noClone){

 using namespace std;

 unsigned int iThread = omp_get_thread_num();
 if(name.find("Template")!=std::string::npos) iThread = 0;

  std::unordered_map<std::string,TH3F*>::iterator it = my3Dhistograms_[iThread].find(name);

 if(noClone && it!=my3Dhistograms_[iThread].end()) it->second;
 if(it!=my3Dhistograms_[iThread].end()){
   TH3F* hClone =  (TH3F*)(it->second->Clone());
   hClone->SetDirectory(0);
   return hClone;
 }
 else cout<<"ERROR: Histogram: "<<name<<" not found in thread: "<<iThread<<endl;
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
void AnalysisHistograms::finalizeHistograms(){

  for(unsigned int iThread = 1;iThread<omp_get_max_threads();++iThread){
    for(auto it:my1Dhistograms_[0]){
      if(my1Dhistograms_[iThread].find(it.first)!=my1Dhistograms_[iThread].end()) it.second->Add(my1Dhistograms_[iThread].find(it.first)->second);
    }
    for(auto it:my2Dhistograms_[0]){
      if(my2Dhistograms_[iThread].find(it.first)!=my2Dhistograms_[iThread].end()) it.second->Add(my2Dhistograms_[iThread].find(it.first)->second);      
    }
    for(auto it:my3Dhistograms_[0]){
      if(my3Dhistograms_[iThread].find(it.first)!=my3Dhistograms_[iThread].end()) it.second->Add(my3Dhistograms_[iThread].find(it.first)->second);      
    }    
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::clear(){

  unsigned int iThread = omp_get_thread_num();
  
  for_each(myProfiles_[iThread].begin(), myProfiles_[iThread].end(), &AnalysisHistograms::resetHistos);  
  for_each(my1Dhistograms_[iThread].begin(), my1Dhistograms_[iThread].end(), &AnalysisHistograms::resetHistos);
  for_each(my2Dhistograms_[iThread].begin(), my2Dhistograms_[iThread].end(), &AnalysisHistograms::resetHistos);
  for_each(my3Dhistograms_[iThread].begin(), my3Dhistograms_[iThread].end(), &AnalysisHistograms::resetHistos);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::write(){

  finalizeHistograms();

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
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
