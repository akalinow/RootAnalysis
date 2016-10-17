#!/usr/bin/env python

from ROOT import *


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV.root"

WAW_file = TFile(WAW_fileName)

SYNCH_file = TFile(SYNCH_fileName,"RECREATE")
inclusive_dir = SYNCH_file.mkdir("mt_inclusive")


#WAW to SYNCH histograms names map
histoPrefix = "HTTAnalyzer/h1DMassVis"
#histoPrefix = "HTTAnalyzer/h1DPtMET"

histogramsMap = {
    histoPrefix+"Data":"data_obs",
    histoPrefix+"DYMuTauJets":"ZTT",
    histoPrefix+"DYMuMuJets":"ZL",
    histoPrefix+"DYOtherJets":"ZJ",
    histoPrefix+"WJets":"W",
    histoPrefix+"TTbar":"TT",
    histoPrefix+"ST":"T",
    histoPrefix+"DiBoson":"VV",
    histoPrefix+"QCDEstimate":"QCD",
    histoPrefix+"ggH120":"ggH120",
    histoPrefix+"qqH120":"qqH120",
    histoPrefix+"ggH125":"ggH125",
    histoPrefix+"qqH125":"qqH125",
    histoPrefix+"ggH130":"ggH130",
    histoPrefix+"qqH130":"qqH130"
    }

for key,value in histogramsMap.items():
    print key        
    histogram = WAW_file.Get(key)                
    histogram.SetName(value)
    histogram.SetDirectory(inclusive_dir)

SYNCH_file.Write()
