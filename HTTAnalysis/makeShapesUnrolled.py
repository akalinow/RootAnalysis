#!/usr/bin/env python

from ROOT import *
import array


WAW_fileName = "RootAnalysis_Analysis.root"
SYNCH_fileName = "htt_mt.inputs-sm-13TeV-2D.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = {
         "mt_0jet":"HTTAnalyzer/h1DUnRollTauPtMassVis",
         "mt_boosted":"HTTAnalyzer/h1DUnRollHiggsPtMassSV",
         "mt_vbf":"HTTAnalyzer/h1DUnRollMjjMassSV"
         }

#define number of bins in each category (nbinsX, nbinsY)
nbins = {
         "mt_0jet":(13,7),
         "mt_boosted":(11,7),
         "mt_vbf":(6,5)
         }

categoryNames = [
                 "muTau_0jet_low", "muTau_0jet_high",
                 "muTau_1jet_low", "muTau_1jet_high",
                 "muTau_vbf_low",  "muTau_vbf_high",
                 "mt_0jet", "mt_0jetCP", 
                 "mt_boosted", "mt_vbf" 
                 ]

relevantCategories = [
                "mt_0jet", 
                 "mt_boosted", "mt_vbf" 
                 ]

histogramsMap = {
    "Data_OS":"data_obs",
    "DYJetsMatchT_OS":"ZTT",
    "DYJetsMatchL_OS":"ZL",
    "DYJetsMatchJ_OS":"ZJ",
    "WJets_OS":"W",
    "TTbarMatchJ_OS":"TTJ",
    "TTbarMatchT_OS":"TTT",
    "ST_OS":"T",
    "DiBosonMatchT_OS":"VVT",
    "DiBosonMatchJ_OS":"VVJ",
    "QCDEstimate":"QCD",
    "ggH120_OS":"ggH120",
    "qqH120_OS":"qqH120",
    "ggH125_OS":"ggH125",
    "qqH125_OS":"qqH125",
    "ggH130_OS":"ggH130",
    "qqH130_OS":"qqH130",
    "ZH120_OS":"ZH120",
    "WplusH120_OS":"WH120",
    "ZH125_OS":"ZH125",
    "WplusH125_OS":"WH125",
    "ZH130_OS":"ZH130",
    "WplusH130_OS":"WH130",
    "EWK2Jets_OS":"EWKZ"
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

channel="mt_"

import array

hData = 0
for iCategory in xrange(0,len(categoryNames)):
    categoryName = categoryNames[iCategory]
    if categoryName not in relevantCategories: continue
    categoryDir = SYNCH_file.mkdir(categoryName)

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            histogram = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
            #continue
        histogram.SetName(value)
        
        histogram.SetDirectory(categoryDir)
        if value=="data_obs": continue
        
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
                    
                    hNameUp = histoPrefix[categoryName] + key+"_"+str(iCategory)+nuisanceParam+"Up"
                    histogramUp = WAW_file.Get(hNameUp)
                    if(histogramUp==None):
                        histogramUp = histogram.Clone()
                        #histogramUp.Reset()
                    histogramUp.SetName(value+"_"+nuisanceParam+"Up")
                    
                    hNameDown = histoPrefix[categoryName] + key+"_"+str(iCategory)+nuisanceParam+"Down"
                    histogramDown = WAW_file.Get(hNameDown)
                    if(histogramDown==None):
                        histogramDown = histogram.Clone()
                        #histogramDown.Reset()
                    histogramDown.SetName(value+"_"+nuisanceParam+"Down")
                    
                    histogramUp.SetDirectory(categoryDir)
                    histogramDown.SetDirectory(categoryDir)
                    SYNCH_file.Write()
                    
                    if cat=="": break
                if (chn=="" or done1): break
        
        #SYNCH_file.Write()

#SYNCH_file.Write()
