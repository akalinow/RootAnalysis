#ifndef commonUtils_H
#define commonUtils_H

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

using namespace std;
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//Return a number in nice format using 10^x powers.
//to be used in making the tables with event counts
string expoEff(float number, float error=0);

////////////////////////////////////////////////////
///Calculate error on the rejection = 1/efficiency
//using binomial formula
double rejError(float rejection, int nEvents);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///Evaluate a value on TGraph, by finding closest
//value on X axis. No interpolation is made.
//if error==true a corresponding error value is returned
double Eval(TGraph* gr, float val, bool error=false);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///Return integrated version of histogram.
///Integration runs from bin 0.
///Method expects TH1F type
TH1F * Integrate(TH1F * histoD);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///Return integrated version of histogram.
///Integration runs from bin 0.
///Method expects TH1D type
TH1D * Integrate(TH1D * histoD);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///Return integrated version of histogram.
///if opt==1 the integration runs diagonally from high end
///if opt==0 the integration runs diagonally from 0
TH2F *Integrate(TH2F *histoD, int opt=0);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//Transform TGraph containing efficiency values to
//rejection=1/eff 
TGraph *grEffToRej(TGraph *grEff);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///Make TGraph containin the ROC curve.
///hSgn and hBkg contain distributions of the
///discriminating variable
///opt=1 returns TGraph with reversed order of points
TGraph* getSgnVsBkg(TH1F *hSgn, TH1F *hBkg, int opt=0);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///Return TCanvas with default setting on margins and background
TCanvas * getDefaultCanvas(float x=10,float y=30,float w=650,float h=600);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///Make a projection of 2D histogramd to 1D.
/// file is a pointer to a TFile containing the 2D histogram
/// proj=X, Y selects axis on which the projection is made
/// low and high parameters control range of projection
TH1F* get1DHisto(string hName,string proj,
		 float low, float high,
		 TFile *file);
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///Set default TLegend properties
void setupLegend(TLegend *leg);

#endif
