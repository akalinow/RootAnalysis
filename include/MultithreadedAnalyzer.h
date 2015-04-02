

#ifndef MULTITHREADEDANALYZER_H
#define MULTITHREADEDANALYZER_H
#include"Analyzer.h"
#include"AnalysisHistograms.h"
#include"../OTFAnalysis/OTFAnalyzer.h"

#include"../OTFAnalysis/OTFHistograms.h"

#include "omp.h"

#include "TH2.h"
#include<vector>
using namespace std;
class MultithreadedAnalyzer : public Analyzer
{
int threads = 16;
public:
    MultithreadedAnalyzer(const string& aName, int threadnum);
    virtual ~MultithreadedAnalyzer();
    
    virtual void initialize(TFileDirectory& aDir, pat::strbitset* aSelections);
    virtual bool analyze(const EventProxyBase& iEvent);
    virtual void finalize();
    virtual void addBranch(TTree* tree);
    virtual void clear();
    virtual void addCutHistos(TList* aList);
private:
    vector<Analyzer*> analyzers;
    vector<AnalysisHistograms*> histograms; 
    AnalysisHistograms * finalHistos;
   void piece1DHistogramsTogether();
   void piece2DHistogramsTogether();
   void piece3DHistogramsTogether();
};

#endif // MULTITHREADEDANALYZER_H
