#ifndef SVfitHistograms_h
#define SVfitHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

#include "HTTEvent.h"
#include "AnalysisHistograms.h"

class THStack;

class SVfitHistograms : public AnalysisHistograms {
public:

SVfitHistograms(TDirectory *myDir);

SVfitHistograms(TDirectory *myDir, const std::vector<std::string> & flavours, std::string channel = "");

virtual ~SVfitHistograms();

void finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester);

using AnalysisHistograms::get1DHistogram;

std::string getTemplateName(const std::string& name);

private:

virtual void defineHistograms();

///Types of the selection flow
std::vector<std::string> selectionFlavours_;

///Plot the ROC curve for given variable(hName) and signal and background processes. 
TH1D* plotROC(std::string hName, std::string signal, std::string background);

///Draw the ROC curves for many discriminators on a single plot and signal and background processes. 
void plotROC(std::string signal, std::string background);
 

//Plot a single 1D histogram. One has to provide the full
//histogram name, e.g. including h1D prefix.
void plotSingleHistogram(std::string hName);

//Plot a single 2D histogram. One has to provide the full
//histogram name, e.g. including h2D prefix.
void plotSingleHistogram2D(std::string hName);


 void plot3DProfile(std::string hName, std::string option = "xy");

std::vector<const HTTAnalysis::eventCategory*> myCategoryRejester;

std::stringstream outputStream;
std::string myChannel_;

};

#endif
