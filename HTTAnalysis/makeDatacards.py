#!/usr/bin/env python

from ROOT import *


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV.root"

WAW_file = TFile(WAW_fileName)

SYNCH_file = TFile(SYNCH_fileName,"RECREATE")
inclusive_dir = SYNCH_file.mkdir("mt_inclusive")


#WAW to SYNCH histograms names map
histoPrefix = "HTTAnalyzer/h1DMassVis"

histogramsMap = {
    histoPrefix+"Data":"data_obs",
    histoPrefix+"DYJetsMuTau":"ZTT",
    histoPrefix+"WJets":"W",
    histoPrefix+"TTbar":"TT",
    histoPrefix+"QCDEstimate":"QCD",
    histoPrefix+"ggH":"ggH125",
    histoPrefix+"VBFH":"vbfH125"
    }

for key,value in histogramsMap.items():
    histogram = WAW_file.Get(key)
    histogram.SetName(value)
    histogram.SetDirectory(inclusive_dir)

SYNCH_file.Write()
