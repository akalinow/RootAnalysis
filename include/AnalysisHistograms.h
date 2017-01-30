#ifndef PFAnalyses_CommonTools_AnalysisHistograms_h
#define PFAnalyses_CommonTools_AnalysisHistograms_h

// Base class for histogram managing.
//
// Original Author:  Artur Kalinowski
//         Created:  Tue Oct 24 15:06:12 CEST 2006
// $Id: AnalysisHistograms.h,v 1.12 2010/09/14 12:20:58 cbern Exp $
//
//
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <vector>
#include <omp.h>

#include "TFile.h"

//ROOT include files
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

class AnalysisHistograms {
 public:

  AnalysisHistograms():file_(0), histosInitialized_(false){};


  AnalysisHistograms(TDirectory *myDir,
		     const std::string & name=""): file_(0),
    histosInitialized_(false), name_(name){};

  static const unsigned int maxThreads = 128;

  virtual ~AnalysisHistograms();

  virtual bool fillProfile(const std::string& name, float x, float val, float weight=1.0);
  virtual bool fill1DHistogram(const std::string &name, float val, float weight=1.0);
  virtual bool fill2DHistogram(const std::string &name, float val1, float val2, float weight=1.0);
  virtual bool fill3DHistogram(const std::string &name, float val1, float val2, float val3, float weight=1.0);

  ///Fill 1D histogram with unrolled Y direction.
  ///Bins from each row in Y are set side by side in 1D.
  virtual bool fill2DUnrolledHistogram(const std::string &name, float val1, float val2, float weight=1.0);

  TProfile* getProfile(const std::string& name, bool noClone = false);
  TH1F* get1DHistogram(const std::string& name, bool noClone = false);
  TH2F* get2DHistogram(const std::string& name, bool noClone = false);
  TH3F* get3DHistogram(const std::string& name, bool noClone = false);

  virtual void finalizeHistograms();

  virtual std::string getTemplateName(const std::string & hName) {return "";};

 protected:

  void clear();

  void write();

  virtual void init(TDirectory *myDir, const std::string & name="");

  virtual void init(const std::string & name=""){};

  virtual void defineHistograms() = 0;

  /// The ROOT file with histograms
  TDirectory *file_;

  ///The secondary directory.
  std::vector<TDirectory*> mySecondaryDirs_;

  /// The histograms
  std::unordered_map<std::string,TProfile*> myProfiles_[maxThreads];
  std::unordered_map<std::string,TH1F*> my1Dhistograms_[maxThreads];
  std::unordered_map<std::string,TH2F*> my2Dhistograms_[maxThreads];
  std::unordered_map<std::string,TH3F*> my3Dhistograms_[maxThreads];

  void addProfile(const std::string& name, const std::string& title,
		  int nBinsX, float xlow, float xhigh,
		  TDirectory* myDir);

  void addRollHistogram(const std::string& name, const std::string& title,
			const std::vector<double> & binsX,
			const std::vector<double> & binsY,
			TDirectory* myDir);

  void add1DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, float xlow, float xhigh,
		      TDirectory* myDir);

  void add1DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, const double* bins,
		      TDirectory* myDir);

  void add2DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, float xlow, float xhigh,
		      int nBinsY, float ylow, float yhigh,
		      TDirectory* myDir);

  void add2DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, const double* binsX,
		      int nBinsY, const double* binsY,
		      TDirectory* myDir);

  void add2DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, float xlow, float xhigh,
		      int nBinsY, const double* binsY,
		      TDirectory* myDir);

  void add3DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, float xlow, float xhigh,
		      int nBinsY, float ylow, float yhigh,
		      int nBinsZ, float zlow, float zhigh,
		      TDirectory* myDir);

  void add3DHistogram(const std::string& name, const std::string& title,
		      int nBinsX, const double* binsX,
		      int nBinsY, const double* binsY,
		      int nBinsZ, const double* binsZ,
		      TDirectory* myDir);

  static void resetHistos(std::pair<const std::string, TH1*> aPair);

  TDirectory *myDirCopy;
  double* equalRanges(int nSteps, double min, double max, double *ranges);
  bool histosInitialized_;

  ///Name of the histograms set instance
  std::string name_;

};

#endif
