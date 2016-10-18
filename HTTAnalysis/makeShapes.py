#!/usr/bin/env python

from ROOT import *


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = "HTTAnalyzer/h1DMassVis"
#histoPrefix = "HTTAnalyzer/h1DPtMET"

categoryNames = ["muTau_0jet_low", "muTau_0jet_high",
                 "muTau_1jet_low", "muTau_1jet_high",
                 "muTau_vbf_low",  "muTau_vbf_high"]

histogramsMap = {
    histoPrefix+"Data_OS":"data_obs",
    histoPrefix+"DYZTTJets_OS":"ZTT",
    histoPrefix+"DYZLJets_OS":"ZL",
    histoPrefix+"DYZJJets_OS":"ZJ",
    histoPrefix+"WJets_OS":"W",
    histoPrefix+"TTbar_OS":"TT",
    histoPrefix+"ST_OS":"T",
    histoPrefix+"DiBoson_OS":"VV",
    histoPrefix+"QCDEstimate":"QCD",
    histoPrefix+"ggH120_OS":"ggH120",
    histoPrefix+"qqH120_OS":"qqH120",
    histoPrefix+"ggH125_OS":"ggH125",
    histoPrefix+"qqH125_OS":"qqH125",
    histoPrefix+"ggH130_OS":"ggH130",
    histoPrefix+"qqH130_OS":"qqH130"
    }


hData = 0
for iCategory in xrange(0,6):
    categoryDir = SYNCH_file.mkdir(categoryNames[iCategory])
    for key,value in histogramsMap.items():
        hName = key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"

        histogram.SetName(value)
        histogram.SetDirectory(categoryDir)

SYNCH_file.Write()
