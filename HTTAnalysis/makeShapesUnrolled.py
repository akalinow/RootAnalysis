#!/usr/bin/env python

from ROOT import *
import array


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV-2D.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = "HTTAnalyzer/h1DUnRollTauPtMassVis"

categoryNames = [
		 "muTau_0jet_low", "muTau_0jet_high",
                 "muTau_1jet_low", "muTau_1jet_high",
                 "muTau_vbf_low",  "muTau_vbf_high",
		 "mt_jet0", "mt_boosted", "mt_vbf" ]

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
    
#according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#Systematic_uncertainties
nuisanceParams = {
    "CMS_shape_t_mt_13TeV",
    "CMS_shape_t_tt_13TeV",
    #"CMS_scale_e_em_13TeV",
    "CMS_scale_j_13TeV",
    "CMS_htt_jetToTauFake_13TeV",
    "CMS_htt_ZLShape_mt_13TeV",
    "CMS_htt_ZLShape_mt_13TeV",
    "CMS_htt_dyShape_13TeV",
    "CMS_htt_ttbarShape_13TeV",
    "QCDSFUncert_mt_0jet_13TeV",
    "QCDSFUncert_mt_boosted_13TeV",
    "QCDSFUncert_mt_vbf_13TeV",
    "WSFUncert_mt_0jet_13TeV",
    "WSFUncert_mt_boosted_13TeV",
    "WSFUncert_mt_vbf_13TeV",
    "CMS_scale_gg_13Tev",
    "CMS_htt_zmumuShape_VBF_13TeV",
    "CMS_htt_zmumuShape_boosted_13TeV"
    }

import array
inclusiveBinning = array.array('d', [0,50,100,300,500,1000])

hData = 0
for iCategory in xrange(6,9):
    categoryName = categoryNames[iCategory]
    categoryDir = SYNCH_file.mkdir(categoryName)

    for key,value in histogramsMap.items():
        hName = key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            continue
        histogram.SetName(value)
            
        for nuisanceParam in nuisanceParams:
        
            hNameUp = key+"_"+str(iCategory)+nuisanceParam+"Up"
            histogramUp = WAW_file.Get(hNameUp)
            if(histogramUp==None):
                histogramUp = histogram.Clone()
                histogramUp.Reset()
            histogramUp.SetName(value+"_"+nuisanceParam+"Up")
            
            hNameDown = key+"_"+str(iCategory)+nuisanceParam+"Down"
            histogramDown = WAW_file.Get(hNameDown)
            if(histogramDown==None):
                histogramDown = histogramUp
            histogramDown.SetName(value+"_"+nuisanceParam+"Down")
            
            histogramUp.SetDirectory(categoryDir)
            histogramDown.SetDirectory(categoryDir)

        #officialBinning = array.array('d', binningMap[categoryName])
        #histogramRebin = histogram.Rebin(len(officialBinning)-1, value, officialBinning)
        #histogramRebin.SetDirectory(categoryDir)
        histogram.SetDirectory(categoryDir)

SYNCH_file.Write()
