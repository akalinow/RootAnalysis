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

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void AnalysisHistograms::addProfile(const std::string& name,
                                    const std::string& title,
                                    int nBinsX, float xlow, float xhigh,
                                    TDirectory* myDir) {

        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TProfile *hTmp = 0;
        hTmp = new TProfile(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);

        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(myProfiles_[iThread].find(name)==myProfiles_[iThread].end()) myProfiles_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing profile!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::addRollHistogram(const std::string& name, const std::string& title,
                                          const std::vector<double> & binsX,
                                          const std::vector<double> & binsY,
                                          TDirectory* myDir){


        int nBinsX = binsX.size()-1;
        int nBinsY = binsY.size()-1;
        ///Include only overflow bins in unrolled histogram.
        int nUnrolledBins = nBinsX*nBinsY;

        TString rolledName(name.c_str());
        rolledName.ReplaceAll("UnRoll","Roll");
        rolledName.ReplaceAll("1D","2D");

        add2DHistogram(rolledName.Data(),title,
                       nBinsX, binsX.data(),
                       nBinsY, binsY.data(),
                       file_);
        add1DHistogram(name,title,nUnrolledBins, 0, nUnrolledBins,file_);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, float xlow, float xhigh,
                                        TDirectory* myDir){

        using namespace std;
        unsigned int iThread = omp_get_thread_num();
        TH1F *hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my1Dhistograms_[iThread].find(name)==my1Dhistograms_[iThread].end()) my1Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add1DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, const double* bins,
                                        TDirectory* myDir){

        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH1F *hTmp = new TH1F(name.c_str(),title.c_str(),nBinsX,bins);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my1Dhistograms_[iThread].find(name)==my1Dhistograms_[iThread].end()) my1Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, float xlow, float xhigh,
                                        int nBinsY, float ylow, float yhigh,
                                        TDirectory* myDir){

        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH2F *hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh, nBinsY, ylow,yhigh);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, const double* binsX,
                                        int nBinsY, const double* binsY,
                                        TDirectory* myDir){
        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH2F *hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add2DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, float xlow, float xhigh,
                                        int nBinsY, const double* binsY,
                                        TDirectory* myDir){

        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH2F *hTmp = new TH2F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,binsY);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my2Dhistograms_[iThread].find(name)==my2Dhistograms_[iThread].end()) my2Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, float xlow, float xhigh,
                                        int nBinsY, float ylow, float yhigh,
                                        int nBinsZ, float zlow, float zhigh,
                                        TDirectory* myDir){
        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH3F *hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,xlow,xhigh,nBinsY,ylow,yhigh, nBinsZ,zlow,zhigh);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        if(my3Dhistograms_[iThread].find(name)==my3Dhistograms_[iThread].end()) my3Dhistograms_[iThread][name] = hTmp;
        else cout<<"ERROR Substituting existing histogram!"<<endl;
}
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::add3DHistogram(const std::string& name, const std::string& title,
                                        int nBinsX, const double* binsX,
                                        int nBinsY, const double* binsY,
                                        int nBinsZ, const double* binsZ,
                                        TDirectory* myDir){

        using namespace std;
        unsigned int iThread = omp_get_thread_num();

        TH3F *hTmp = new TH3F(name.c_str(),title.c_str(),nBinsX,binsX,nBinsY,binsY,nBinsZ,binsZ);

        hTmp->Sumw2();
        hTmp->SetDirectory(0);

        //TEST if(iThread==0 && name.find("Template")==std::string::npos) hTmp->SetDirectory(myDir);
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

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
        else{
                std::string templateName = getTemplateName(name);
                TProfile *pTemplate = getProfile(templateName,true);
                if(!pTemplate) return false;

                addProfile(name,"",
                           pTemplate->GetNbinsX(),
                           pTemplate->GetXaxis()->GetXmin(),
                           pTemplate->GetXaxis()->GetXmax(),
                           file_);

                myProfiles_[iThread].find(name)->second->Fill(x,val,weight);
        }
        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill1DHistogram(const std::string& name, float val, float weight){

        using namespace std;

        unsigned int iThread = omp_get_thread_num();

        std::unordered_map<std::string,TH1F*>::iterator it = my1Dhistograms_[iThread].find(name);

        if(it!=my1Dhistograms_[iThread].end()) it->second->Fill(val,weight);
        else{
                std::string templateName = getTemplateName(name);
                TH1F *hTemplate = get1DHistogram(templateName,true);
                if(!hTemplate) return false;
                if(hTemplate->GetXaxis()->IsVariableBinSize()) {
                        add1DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXbins()->GetArray(),
                                       file_);
                }
                else{
                        add1DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXmin(),
                                       hTemplate->GetXaxis()->GetXmax(),
                                       file_);
                }
                my1Dhistograms_[iThread].find(name)->second->Fill(val,weight);
        }
        return true;
}
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

        using namespace std;

        unsigned int iThread = omp_get_thread_num();

        std::unordered_map<std::string,TH2F*>::iterator it = my2Dhistograms_[iThread].find(name);
        if(it!=my2Dhistograms_[iThread].end()) it->second->Fill(val1,val2,weight);
        else{
                std::string templateName = getTemplateName(name);
                TH2F *hTemplate = get2DHistogram(templateName,true);
                if(!hTemplate) return false;
                if(hTemplate->GetXaxis()->IsVariableBinSize() &&
                   hTemplate->GetYaxis()->IsVariableBinSize()) {
                        add2DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXbins()->GetArray(),
                                       hTemplate->GetNbinsY(),
                                       hTemplate->GetYaxis()->GetXbins()->GetArray(),
                                       file_);
                }
                else{
                        add2DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXmin(),
                                       hTemplate->GetXaxis()->GetXmax(),
                                       hTemplate->GetNbinsY(),
                                       hTemplate->GetYaxis()->GetXmin(),
                                       hTemplate->GetYaxis()->GetXmax(),
                                       file_);
                }
                my2Dhistograms_[iThread].find(name)->second->Fill(val1,val2,weight);
        }

        return true;
}
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill3DHistogram(const std::string& name, float val1, float val2, float val3,
                                         float weight){

        using namespace std;

        unsigned int iThread = omp_get_thread_num();

        std::unordered_map<std::string,TH3F*>::iterator it = my3Dhistograms_[iThread].find(name);
        if(it!=my3Dhistograms_[iThread].end()) it->second->Fill(val1,val2,val3,weight);
        else{
                std::string templateName = getTemplateName(name);
                TH3F *hTemplate = get3DHistogram(templateName,true);
                if(!hTemplate) return false;
                if(hTemplate->GetXaxis()->IsVariableBinSize()) {
                        add3DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXbins()->GetArray(),
                                       hTemplate->GetNbinsY(),
                                       hTemplate->GetYaxis()->GetXbins()->GetArray(),
                                       hTemplate->GetNbinsZ(),
                                       hTemplate->GetZaxis()->GetXbins()->GetArray(),
                                       file_);
                }
                else{
                        add3DHistogram(name,"",
                                       hTemplate->GetNbinsX(),
                                       hTemplate->GetXaxis()->GetXmin(),
                                       hTemplate->GetXaxis()->GetXmax(),
                                       hTemplate->GetNbinsY(),
                                       hTemplate->GetYaxis()->GetXmin(),
                                       hTemplate->GetYaxis()->GetXmax(),
                                       hTemplate->GetNbinsZ(),
                                       hTemplate->GetZaxis()->GetXmin(),
                                       hTemplate->GetZaxis()->GetXmax(),
                                       file_);
                }
                my3Dhistograms_[iThread].find(name)->second->Fill(val1,val2,val3,weight);
        }

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool AnalysisHistograms::fill2DUnrolledHistogram(const std::string &name, float val1, float val2, float weight){

        bool noClone = true;
        TString rolledName(name.c_str());
        rolledName.ReplaceAll("UnRoll","Roll");
        rolledName.ReplaceAll("1D","2D");
        fill2DHistogram(rolledName.Data(), val1, val2, weight);
        TH2F *hRolled = get2DHistogram(rolledName.Data(), noClone);

        if(!hRolled) return false;
        int nBinsX = hRolled->GetNbinsX();
        int iUnrolledBinX = hRolled->GetXaxis()->FindBin(val1);
        int iUnrolledBinY = hRolled->GetYaxis()->FindBin(val2);
        int iUnrolledBin = iUnrolledBinX + (iUnrolledBinY-1)*nBinsX;

        return fill1DHistogram(name, iUnrolledBin-0.5, weight);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TProfile* AnalysisHistograms::getProfile(const std::string& name, bool noClone) {

        using namespace std;

        unsigned int iThread = omp_get_thread_num();
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        std::unordered_map<std::string,TProfile*>::iterator it = myProfiles_[iThread].find(name);

        if(noClone && it!=myProfiles_[iThread].end()) return it->second;
        else if(it!=myProfiles_[iThread].end()) {
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
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        std::unordered_map<std::string,TH1F*>::iterator it = my1Dhistograms_[iThread].find(name);

        if(noClone && it!=my1Dhistograms_[iThread].end()) return it->second;
        else if(it!=my1Dhistograms_[iThread].end()) {
                TH1F* hClone =  (TH1F*)(it->second->Clone());
                hClone->SetDirectory(0);
                return hClone;
        }
        else {
                ///To many ERRORS for histograms absent due to samples not passing the selection.
                //cout<<"ERROR: Histogram: "<<name<<" not found in thread: "<<omp_get_thread_num()<<endl;
        }

        return 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TH2F* AnalysisHistograms::get2DHistogram(const std::string& name, bool noClone){

        using namespace std;

        unsigned int iThread = omp_get_thread_num();
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        std::unordered_map<std::string,TH2F*>::iterator it = my2Dhistograms_[iThread].find(name);

        if(noClone && it!=my2Dhistograms_[iThread].end()) return it->second;
        else if(it!=my2Dhistograms_[iThread].end()) {
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
        if(name.find("Template")!=std::string::npos) iThread = AnalysisHistograms::maxThreads-1;

        std::unordered_map<std::string,TH3F*>::iterator it = my3Dhistograms_[iThread].find(name);

        if(noClone && it!=my3Dhistograms_[iThread].end()) it->second;
        if(it!=my3Dhistograms_[iThread].end()) {
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
        if(aPair.second) aPair.second->Reset();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::init(TDirectory *myDir,
                              const std::string & name){

        name_ = name;
        myDirCopy = myDir;
        if(!histosInitialized_) {
                if(name_.size()) {
                        mySecondaryDirs_.push_back(myDir->mkdir(name_.c_str()));
                        file_ = mySecondaryDirs_[mySecondaryDirs_.size()-1];
                }
                else file_ = myDir;
                defineHistograms();
        }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void AnalysisHistograms::finalizeHistograms(){

///Add histograms missing in thread 0 map.
        for(int iThread = 1; iThread<omp_get_max_threads(); ++iThread) {
                for(auto it:my1Dhistograms_[iThread]) {
                        if(my1Dhistograms_[0].find(it.first)==my1Dhistograms_[0].end()) {
                                TH1F *hEmpty = (TH1F*)it.second->Clone();
                                hEmpty->Reset();
                                my1Dhistograms_[0][it.first] = hEmpty;
                        }
                }
                for(auto it:my2Dhistograms_[iThread]) {
                        if(my2Dhistograms_[0].find(it.first)==my2Dhistograms_[0].end()) {
                                TH2F *hEmpty = (TH2F*)it.second->Clone();
                                hEmpty->Reset();
                                my2Dhistograms_[0][it.first] = hEmpty;
                        }
                }
                for(auto it:my3Dhistograms_[iThread]) {
                        if(my3Dhistograms_[0].find(it.first)==my3Dhistograms_[0].end()) {
                                TH3F *hEmpty = (TH3F*)it.second->Clone();
                                hEmpty->Reset();
                                my3Dhistograms_[0][it.first] = hEmpty;
                        }
                }
                for(auto it:myProfiles_[iThread]) {
                        if(myProfiles_[0].find(it.first)==myProfiles_[0].end()) {
                                TProfile *hEmpty = (TProfile*)it.second->Clone();
                                hEmpty->Reset();
                                myProfiles_[0][it.first] = hEmpty;
                        }
                }
        }

        std::cout<<"1D histogram size: "<<my1Dhistograms_[0].size()<<std::endl;
        std::cout<<"2D histogram size: "<<my2Dhistograms_[0].size()<<std::endl;
        std::cout<<"3D histogram size: "<<my3Dhistograms_[0].size()<<std::endl;
        std::cout<<"TProfile size: "<<myProfiles_[0].size()<<std::endl;

        for(int iThread = 1; iThread<omp_get_max_threads(); ++iThread) {
                if(my1Dhistograms_[0].size()!=my1Dhistograms_[iThread].size()) {
                        std::cout<<"1D histogram size mismatch. "
                                 <<" thread 0: "<<my1Dhistograms_[0].size()
                                 <<" thread "<<iThread<<" "<<my1Dhistograms_[iThread].size()
                                 <<std::endl;
                }
                if(my2Dhistograms_[0].size()!=my2Dhistograms_[iThread].size()) {
                        std::cout<<"2D histogram size mismatch. "
                                 <<" thread 0: "<<my2Dhistograms_[0].size()
                                 <<" thread "<<iThread<<" "<<my2Dhistograms_[iThread].size()
                                 <<std::endl;
                }
                if(my3Dhistograms_[0].size()!=my3Dhistograms_[iThread].size()) {
                        std::cout<<"3D histogram size mismatch. "
                                 <<" thread 0: "<<my3Dhistograms_[0].size()
                                 <<" thread "<<iThread<<" "<<my3Dhistograms_[iThread].size()
                                 <<std::endl;
                }
                if(myProfiles_[0].size()!=myProfiles_[iThread].size()) {
                        std::cout<<"TProfile size mismatch. "
                                 <<" thread 0: "<<myProfiles_[0].size()
                                 <<" thread "<<iThread<<" "<<myProfiles_[iThread].size()
                                 <<std::endl;
                }
                for(auto it:my1Dhistograms_[0]) {
                        auto itThread = my1Dhistograms_[iThread].find(it.first);
                        if(itThread!=my1Dhistograms_[iThread].end()) {
                                it.second->Add(itThread->second);
                                delete itThread->second;
                                my1Dhistograms_[iThread].erase(itThread);
                        }
                }

                for(auto it:my2Dhistograms_[0]) {
                        auto itThread = my2Dhistograms_[iThread].find(it.first);
                        if(itThread!=my2Dhistograms_[iThread].end()) {
                                it.second->Add(itThread->second);
                                delete itThread->second;
                                my2Dhistograms_[iThread].erase(itThread);
                        }
                }
                for(auto it:my3Dhistograms_[0]) {
                        auto itThread = my3Dhistograms_[iThread].find(it.first);
                        if(itThread!=my3Dhistograms_[iThread].end()) {
                                it.second->Add(itThread->second);
                                delete itThread->second;
                                my3Dhistograms_[iThread].erase(itThread);
                        }
                }
                for(auto it:myProfiles_[0]) {
                        auto itThread = myProfiles_[iThread].find(it.first);
                        if(itThread!=myProfiles_[iThread].end()) {
                                it.second->Add(itThread->second);
                                delete itThread->second;
                                myProfiles_[iThread].erase(itThread);
                        }
                }
        }
/*
        for(auto it:my1Dhistograms_[0]) delete it.second;
        for(auto it:my2Dhistograms_[0]) delete it.second;
        for(auto it:myProfiles_[0]) delete it.second;
                std::cout<<"After delete"<<std::endl;
                std::cin.ignore();
*/
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
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double* AnalysisHistograms::equalRanges(int nSteps, double min, double max, double *ranges){

        for(int i=0; i<=nSteps; ++i) {
                double val = min + i*(max-min)/nSteps;
                ranges[i] = val;
        }

        return ranges;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
