#ifndef HTTHistograms_h
#define HTTHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

#include "HTTEvent.h"
#include "AnalysisHistograms.h"

class THStack;

class HTTHistograms : public AnalysisHistograms {
public:

HTTHistograms(TDirectory *myDir);

HTTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours, std::string channel = "");

virtual ~HTTHistograms();

void finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester);

using AnalysisHistograms::get1DHistogram;

TH1F *get1DHistogram(unsigned int iCategory, std::string varName,
                     unsigned int iSystEffect = (unsigned int) HTTAnalysis::NOMINAL);

std::string getTemplateName(const std::string& name);

///Return sample normalisation based only on
///luminosity and cross section.
///MC to DATA scaling factors should be applied
///on top of this normalisation.
float getSampleNormalisation(std::string sampleName);

///Return luminosity corresponding to given MC sample.
float getSampleLuminosity(const std::string& sampleName, float crossSection);

///Estimate QCD background using the SS/OS method.
TH1F* getQCDbackground(unsigned int iCategory, std::string varName,
                       unsigned int iSystEffect = (unsigned int) HTTAnalysis::NOMINAL);

///Calculate scaling factor for the WJets MC
///Scaling factor is estimated in high Mt region.
///Other backgrounds are subtracted, basing on MC
///QCD contribution is neglected.
std::pair<float,float> getWNormalisation(unsigned int iCategory,
                                         unsigned int iSystEffect = (unsigned int) HTTAnalysis::NOMINAL);

///Calculate QCD ratio between signal and control regions.
std::pair<float,float> getQCDControlToSignal(unsigned int iCategory,
                                             unsigned int iSystEffect = (unsigned int) HTTAnalysis::NOMINAL);

private:

virtual void defineHistograms();

///Types of the selection flow
std::vector<std::string> selectionFlavours_;

//Plot stacked histograms for each contributing process.
///iCategory - selection category to be plotted.
///iSystEffect - systematic effect type to be plotted
///selName - secondary type of selection (OS/SS/mt) used for background estimation
//varName - name of variable to be plotted,
THStack* plotStack(unsigned int iCategory,
                   std::string varName,
                   unsigned int iSystEffect = (unsigned int) HTTAnalysis::NOMINAL);

void plotnPCA(const std::string & type);

void plotVerticesPulls(const std::string & hName);

void plotProfiles(const std::string & hName,
                  const std::string & sysType);

void plotPhiDecayPlanes(const std::string& name);

void plot_HAZ_Histograms(const std::string & hName,
                         const std::string & sysType);

void plotCPhistograms(unsigned int iCategory);

///Return sum of all non Higgs MC contributions.
///sumForW control W MC to DATA correction. By default
///the correction is estimated. To avoid circular dependency,
///it has to be switched off when calculating the scale factor itself.
///as getWNormalisation() method calls this one.
TH1F *getMCSum(unsigned int iCategory, std::string varName, unsigned int iSystEffect, bool sumForW = false);

///Return histogram for sum of all DY decay modes, and jet bins
TH1F *get1D_DYJet_Histogram(const std::string& name);

///Return histogram for sum of all jet bins
TH1F *get1D_WJet_Histogram(const std::string& name);

///Return histogram for sum VV processes
TH1F *get1D_VV_Histogram(const std::string& name, std::string tauMatchSuffix = "");

///Return histogram for sum TT processes
TH1F *get1D_TT_Histogram(const std::string& name, std::string tauMatchSuffix = "");

///Return histogram for sum single top processes
TH1F *get1D_ST_Histogram(const std::string& name);

///Return sum of DY histograms. Sum can run over
///decay modes, jet bins, or both.
TH1F *get1D_TauMatchJetSum(const std::string& name, bool sumDecayModes, bool sumJetBins = false);

//Return sum of histograms where you substitute pattern with all sample names.
//If tauMatchSuffix is specified, you will have sum of histos that match the suffix
TH1F* get1D_SumPattern_Histogram(const std::string& name, const std::string & pattern,
                                 const std::vector<std::string> & sampleNames,
                                 std::string tauMatchSuffix = "");

///Return histogram for sum of all W/Z jet bins
///The results is scaled to 1/LO_xsection.
TH1F *get1D_VJetSum(const std::string& name);

///Return histogram for sum of the EWK + 2jets samples
TH1F *get1D_EWK2JetsSum(const std::string& name);

//Plot a single histogram. One has to provide the full
//histogram name, e.g. including h1D prefix.
void plotSingleHistogram(std::string hName);

std::vector<const HTTAnalysis::eventCategory*> myCategoryRejester;

std::stringstream outputStream;
std::string myChannel_;

};

#endif
