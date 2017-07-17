#include "commonUtils.h"

#include <string>
#include <cmath>
#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"

////////////////////////////////////////////////////
////////////////////////////////////////////////////

using namespace std;

////////////////////////////////////////////////////
////////////////////////////////////////////////////
int colors[14] = {1,2,3,4,8,9,28,29,30,33,38,40,41,42};
////////////////////////////////////////////////////
////////////////////////////////////////////////////
string expoEff(float number, float error){

  char text[200];
  int base;
  if(number>1) base = (unsigned int)(fabs(log(number)/log(10.0)) + 0.0);
  else base = (unsigned int)(fabs(log(number)/log(10.0)) + 0.99);
  if(base<4){
    if(error<0.001) sprintf(text," ${\\mathrm %4.4f \\pm %4.4f}$ ",number,error);
    else if(error<0.01) sprintf(text," ${\\mathrm %3.3f \\pm %3.3f}$ ",number,error);
    else if(error<0.1) sprintf(text," ${\\mathrm %3.2f \\pm %3.2f}$ ",number,error);
    else if(error<10) sprintf(text," ${\\mathrm %3.1f \\pm %3.1f}$ ",number,error);
    else if(error<100) sprintf(text," ${\\mathrm %3.0f \\pm %3.0f}$ ",number,error);
    else if(error<1000) sprintf(text," ${\\mathrm %2.0f \\pm %2.0f}$ ",number,error);
    else if(error<10000) sprintf(text," ${\\mathrm %2.0f \\pm %2.0f}$ ",number,error);
    else std::cout<<"expoEff: error: "<<error<<std::endl;
    string result;
    result.append(text);
    return result;
  }
  int sgn = (unsigned int)(log(number)/fabs(log(number)));
  float tmp = number*pow(10.0,-sgn*base);
  if(base!=0) sprintf(text,"${\\mathrm (%3.2f \\cdot 10^{%d} \\pm",tmp,sgn*base);
  else sprintf(text,"${\\mathrm  %3.2f $",tmp);
  string result;
  result.append(text);

  sgn = (unsigned int)(log(error)/fabs(log(error)));
  tmp = error*pow(10.0,-sgn*base);
  if(base!=0) sprintf(text," %3.2f) \\cdot 10^{%d}} $",tmp,sgn*base);
  else sprintf(text," %3.2f} $",tmp);
  result.append(text);

  return result;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
double rejError(float rejection, int nEvents){

  float eff = 1.0/(rejection+1);
  float effError = sqrt(eff*(1-eff)/nEvents);
  float rejError = effError/eff/eff;
  return rejError;

}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
double Eval(TGraph* gr, float val, bool error){

  int nPoints = gr->GetN();
  Double_t *x;
  Double_t *y;
  x = gr->GetX();
  y = gr->GetY();

  //if(val<x[0]) return gr->Eval(val);
  //else return -999.0;

  double tmp = 0.08;
  int index = -10;
  for(int i=0;i<nPoints;++i){
    if(fabs(x[i]-val)<tmp){
      tmp = fabs(x[i]-val);
      index = i;
    }
  }

  //cout<<"val: "<<val<<" index: "<<index;
  //if(index>0) cout<<" closest: "<<x[index]<<std::endl;

  if(index==-10) return -999.0;
  if(index==0) return -999.0;
  //if(index==0 && tmp>0.01) return -999.0;
  if(error) return gr->GetErrorY(index);
  return y[index];

}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TH1F * Integrate(TH1F * histoD) {

   TH1F * histoI = new TH1F(*histoD);
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
// for (i = 0; i<=histoI->GetNbinsX()+1;i++){
//      cout <<"bin: "<<i<<" cont: "<<histoI->GetBinContent(i);
//      cout <<            " error: "<<histoI->GetBinError(i)<<endl;
// }

   return histoI;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TH1D * Integrate(TH1D * histoD) {

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
// for (i = 0; i<=histoI->GetNbinsX()+1;i++){
//      cout <<"bin: "<<i<<" cont: "<<histoI->GetBinContent(i);
//      cout <<            " error: "<<histoI->GetBinError(i)<<endl;
// }

   return histoI;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TH2F *Integrate(TH2F *histoD, int opt){

  TH2F * histoI = new TH2F(*histoD);
  histoI->Reset();

  float sum = 0.0;
  for (int i = histoD->GetNbinsX()+1; i >= 0; i--) {
    for (int j = histoD->GetNbinsY()+1; j >= 0; j--) {
      if(opt==0)sum = histoD->Integral(i, histoD->GetNbinsX()+1,
				       j,histoD->GetNbinsY()+1);
      if(opt==1 )sum = histoD->Integral(0,i,0,j);
      histoI->SetBinContent(i,j,sum);
    }
  }
  return histoI;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TGraph *grEffToRej(TGraph *grEff){

  int n = grEff->GetN();
  float xx[1000];
  float yy[1000];
  double x,y;

  for(int i=0;i<n;++i){
    grEff->GetPoint(i,x,y);
    xx[i] = x;
    if(y!=0) yy[i] = (1 - y)/y;
    else yy[i] = 100000;
  }

  TGraph *aGraph = new TGraph(n,xx,yy);


  aGraph->SetMarkerStyle(grEff->GetMarkerStyle());
  aGraph->SetMarkerColor(grEff->GetMarkerColor());
  aGraph->SetLineColor(grEff->GetLineColor());

  return aGraph;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TGraph* getSgnVsBkg(TH1F *hSgn, TH1F *hBkg, int opt){

  float x[10000], ex[10000];
  float y[10000], ey[10000];
  float y1[10000], ey1[10000];

  float x2[10000], y2[10000];
  float ex2[10000], ey2[10000];

  TH1F *hSgnInt =  Integrate(hSgn);
  TH1F *hBkgInt =  Integrate(hBkg);

  float lastVal = -1.0;

  int nBins = hSgnInt->GetNbinsX();
  int n = 0;

  for(int i=0;i<=nBins;++i){
    float valSgn = hSgnInt->GetBinContent(i);
    float valBkg = hBkgInt->GetBinContent(i);
    float errSgn = hSgnInt->GetBinError(i);
    float errBkg = hBkgInt->GetBinError(i);
    /////////////////////
    //if(valSgn<0.1) continue;
    //if(nBins>15 && (7*i)%10!=0) continue;
    //if(valBkg==lastVal) continue;
    //if(valBkg!=0) lastVal = valBkg;
    ////////////////////
    x[n] = 1 - valSgn;
    ex[n] = errSgn;
    if(valBkg!=0) {
      y[n] = (1 - valBkg)/valBkg;
      ey[n] =  errBkg/valBkg/valBkg;
    }
    else {
      //continue;
      y[n] = 100000;
      ey[n] = 0;
    }

    y1[n] = 1 - valBkg;
    ey1[n] = errBkg;
    n++;
  }

 for(int i=0;i<n;i++){
    x2[i] = x[n-i-1];
    y2[i] = y1[n-i-1];
    ///////////////////
    ex2[i] = ex[n-i-1];
    ey2[i] = ey1[n-i-1];
  }

 TGraph *aGraph = new TGraph(n,x,y1);
 //TGraph *aGraph1 = new TGraph(n,x2,y2);
 ////Graphs with errors
 //TGraphErrors *aGraph = new TGraphErrors(n,x,y,ex,ey);
 TGraphErrors *aGraph1 = new TGraphErrors(n,x2,y2,ex2,ey2);

  if(opt==1) return aGraph;
  return aGraph1;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
TCanvas * getDefaultCanvas(float x,float y,float w,float h){

  TCanvas *c1 = new TCanvas("c1","Jet LL",x,y,w,h);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetLeftMargin(0.16);
  c1->SetGrid(1,1);
  c1->SetTicky();
  //c1->SetTickx();

  return c1;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void setupLegend(TLegend *leg){

  //leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* get1DHisto(string hName,string proj,
		 float low, float high,
		 TFile *file){

  TH2F *h2D = (TH2F*)file->Get(hName.c_str());

  int binLow = h2D->GetYaxis()->FindBin(low);
  int binHigh = h2D->GetYaxis()->FindBin(high);

  TH1D* h1D = 0;

  std::cout<<"h2D: "<<h2D<<" proj: "<<proj<<std::endl;

  if(proj=="X") h1D = h2D->ProjectionX("h1D",binLow,binHigh);
  if(proj=="Y"){
    binLow = h2D->GetXaxis()->FindBin(low);
    binHigh = h2D->GetXaxis()->FindBin(high);
    h1D = h2D->ProjectionY("h1D",binLow,binHigh);
  }
  if(h1D){
    TH1F *h1DFloat = new TH1F(h1D->GetName(),h1D->GetTitle(),h1D->GetNbinsX(),
			      h1D->GetXaxis()->GetXmin(),
			      h1D->GetXaxis()->GetXmax());
    for(int i=0;i<=h1D->GetNbinsX()+1;++i){
      h1DFloat->SetBinContent(i,h1D->GetBinContent(i));
      h1DFloat->SetBinError(i,h1D->GetBinError(i));
    }
    h1DFloat->Print();
    h1DFloat->SetName("h1F");
    return h1DFloat;
  }

  return 0;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
