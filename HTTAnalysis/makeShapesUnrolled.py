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
                 "mt_0jet", "mt_boosted", "mt_vbf" ]

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
    histoPrefix+"qqH130_OS":"qqH130",
    histoPrefix+"ZH120_OS":"ZH120",
    histoPrefix+"WH120_OS":"WH120",
    histoPrefix+"ZH125_OS":"ZH125",
    histoPrefix+"WH125_OS":"WH125",
    histoPrefix+"ZH130_OS":"ZH130",
    histoPrefix+"WH130_OS":"WH130",
    histoPrefix+"EWK2Jets":"EWKZ"
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
    "CMS_shape_t_13TeV":(("mt_","tt_","et_","em_"),("","")),
    "CMS_scale_t_13TeV":(("mt_","tt_","et_","em_"),("","")),
    "CMS_scale_e_13TeV":(("em_",""),("","")),
    "CMS_scale_j_13TeV":(("",""),("","")),
    "CMS_htt_jetToTauFake_13TeV":(("",""),("","")),
    "CMS_htt_ZLShape_13TeV":(("mt_",""),("","")),
    "CMS_htt_dyShape_13TeV":(("",""),("","")),
    "CMS_htt_ttbarShape_13TeV":(("",""),("","")),
    "QCDSFUncert_13TeV":(("mt_",""),("0jet_","vbf_","boosted_")),
    "WSFUncert_13TeV":(("mt_",""),("0jet_","vbf_","boosted_")),
    "CMS_scale_gg_13TeV":(("",""),("","")),
    "CMS_htt_zmumuShape_13TeV":(("",""),("vbf_","boosted_","VBF_")),
    }

'''
nuisanceParams = {
    "CMS_shape_t_mt_13TeV",
    "CMS_shape_t_tt_13TeV",
    "CMS_scale_j_13TeV",
    "CMS_htt_jetToTauFake_13TeV",
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
'''

channel="mt_"

import array
inclusiveBinning = array.array('d', [0,50,100,300,500,1000])

hData = 0
for iCategory in xrange(6,9):
    categoryName = categoryNames[iCategory]
    categoryDir = SYNCH_file.mkdir(categoryName)

    for key,value in histogramsMap.iteritems():
        hName = key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            histogram = TH1F(value,"",91,0.5,91.5)
            #continue
        histogram.SetName(value)
        
        histogram.SetDirectory(categoryDir)
        
        for key1, value1 in nuisanceParams.items():
            index = key1.find("13TeV")
            name = key1
            for chn in value1[0]:
                if (chn==channel or chn==""): 
                  name = key1[:index]+chn+key1[index:]
                  done1 = True
                  for cat in value1[1]:
                    index = name.find("13TeV")
                    nuisanceParam = name[:index]+cat+name[index:]
                    #print categoryNames[iCategory]+": " + value+"_"+nuisanceParam
                    
                    hNameUp = key+"_"+str(iCategory)+nuisanceParam+"Up"
                    histogramUp = WAW_file.Get(hNameUp)
                    if(histogramUp==None):
                        histogramUp = histogram.Clone()
                        histogramUp.Reset()
                    histogramUp.SetName(value+"_"+nuisanceParam+"Up")
                    
                    hNameDown = key+"_"+str(iCategory)+nuisanceParam+"Down"
                    histogramDown = WAW_file.Get(hNameDown)
                    if(histogramDown==None):
                        histogramDown = histogram.Clone()
                        histogramDown.Reset()
                    histogramDown.SetName(value+"_"+nuisanceParam+"Down")
                    
                    histogramUp.SetDirectory(categoryDir)
                    histogramDown.SetDirectory(categoryDir)
                    SYNCH_file.Write()
                    
                    if cat=="": break
                if (chn=="" or done1): break
        
        #officialBinning = array.array('d', binningMap[categoryName])
        #histogramRebin = histogram.Rebin(len(officialBinning)-1, value, officialBinning)
        #histogramRebin.SetDirectory(categoryDir)
        #SYNCH_file.Write()

#SYNCH_file.Write()
