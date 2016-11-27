#!/usr/bin/env python

from ROOT import *
import array


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = "HTTAnalyzer/h1DMassVis"

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


binningMap = {
    "muTau_0jet_low": range(0,310,10),
    "muTau_0jet_high": range(0,310,10),
    "muTau_1jet_low":  [0,40,60,70,80,90,100,110,120,130,150,200,250],
    "muTau_1jet_high": [0,40,60,70,80,90,100,110,120,130,150,200,250],
    "muTau_vbf_low":   [0,40,60,80,100,120,150,200,250],
    "muTau_vbf_high":  [0,40,60,80,100,120,150,200,250]
    }

import array
inclusiveBinning = array.array('d', [0,50,100,300,500,1000])

hData = 0
for iCategory in xrange(0,6):
    categoryName = categoryNames[iCategory]
    categoryDir = SYNCH_file.mkdir(categoryName)

    for key,value in histogramsMap.items():
        hName = key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            continue

        officialBinning = array.array('d', binningMap[categoryName])
        histogramRebin = histogram.Rebin(len(officialBinning)-1, value, officialBinning)
        histogramRebin.SetDirectory(categoryDir)

SYNCH_file.Write()
