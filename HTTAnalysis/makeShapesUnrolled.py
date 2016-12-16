#!/usr/bin/env python

from ROOT import *
import array


WAW_fileName = "RootAnalysis_Analysis.root"

channel="mt_"

SYNCH_fileName = "htt_mt.inputs-sm-13TeV-2D.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = {
         "mt_0jet":"HTTAnalyzer/h1DUnRollTauPtMassVis",
         "mt_boosted":"HTTAnalyzer/h1DUnRollHiggsPtMassSV",
         "mt_vbf":"HTTAnalyzer/h1DUnRollMjjMassSV"
         }

#define number of bins in each category (nbinsX, nbinsY) in case when you need to create an empty histo
nbins = {
         "mt_0jet":(13,7),
         "mt_boosted":(11,7),
         "mt_vbf":(6,5)
         }

categoryNames = [
                 "muTau_0jet_low", "muTau_0jet_high",
                 "muTau_1jet_low", "muTau_1jet_high",
                 "muTau_vbf_low",  "muTau_vbf_high",

                 "mt_0jet", 
                 "mt_boosted", "mt_vbf",
                 "mt_CP_Phi", "mt_CP_Rho",
                 "mt_wjets_0jet_cr", 
                 "mt_wjets_boosted_cr", "mt_wjets_vbf_cr",
                 "mt_antiiso_0jet_cr",
                 "mt_antiiso_boosted_cr", "mt_antiiso_vbf_cr"
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
    "WH120_OS":"WH120",
    "ZH125_OS":"ZH125",
    "WH125_OS":"WH125",
    "ZH130_OS":"ZH130",
    "WH130_OS":"WH130",
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
    "CMS_htt_zmumuShape_13TeV":(("",""),("0jet_","vbf_","boosted_")),
    }

def getSingleNPHistos(prefix, np, histo):

    hName = prefix+"_"+np
    hUp = WAW_file.Get(hName+"Up")

    if hUp==None :
        print "Missing histogram: ",hName+"Up"
        hUp = histo.Clone()

    hDown = WAW_file.Get(hName+"Down")
    if hDown==None :
        #almost always where there is no Up histo, there is no down histo, so there is no need to print its name again
        #print prefix+np
        hDown = histo.Clone()

    return (hUp, hDown)



def getDoubleNPHistos(prefix, np1, np2, histo):

    hUpUp = WAW_file.Get(prefix+np1+"Up"+np2+"Up")
    if hUpUp==None :
        print "Missing histogram: ",prefix+np1+"Up/Down"+np2+"Up/Down"
        hUpUp = histo.Clone()

    hUpDown = WAW_file.Get(prefix+np1+"Up"+np2+"Down")
    if hUpDown==None :
        hUpDown = histo.Clone()

    hDownUp = WAW_file.Get(prefix+np1+"Down"+np2+"Up")
    if hDownUp==None :
        hDownUp = histo.Clone()

    hDownDown = WAW_file.Get(prefix+np1+"Down"+np2+"Down")
    if hDownDown==None :
        hDownDown = histo.Clone()

    return (hUpUp, hUpDown, hDownUp, hDownDown)

import array

#basic categories
categoryDirMade = False

hData = 0

print "MAIN REGION\n\n\n\n\n"

for iCategory in xrange(0,len(categoryNames)):
    categoryName = categoryNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key+"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            histogram = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram.SetName(value)
        histogram.Write()

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

                    histos = getSingleNPHistos(histoPrefix[categoryName] + key+"_"+str(iCategory), nuisanceParam, histogram)
                    histogramUp = histos[0]
                    histogramUp.SetName(value+"_"+nuisanceParam+"Up")
                    histogramUp.Write()

                    histogramDown = histos[1]
                    histogramDown.SetName(value+"_"+nuisanceParam+"Down")
                    histogramDown.Write()

                    if cat=="": break
                if (chn=="" or done1): break

print "\n\n\n\n\nCONTROL REGIONS\n\n\n\n\n"

##################################################
#control regions
##################################################

histoPrefix = {
         "mt_wjets_0jet_cr":"HTTAnalyzer/h1DMassTrans",
         "mt_wjets_boosted_cr":"HTTAnalyzer/h1DMassTrans",
         "mt_wjets_vbf_cr":"HTTAnalyzer/h1DMassTrans",
         "mt_antiiso_0jet_cr":"HTTAnalyzer/h1DMassVis",
         "mt_antiiso_boosted_cr":"HTTAnalyzer/h1DMassSV",
         "mt_antiiso_vbf_cr":"HTTAnalyzer/h1DMassSV"
         }

histogramsMap = {
    "Data_OS":"data_obs",
    "DYJetsMatchT_OS":"ZTT",
    "DYJetsMatchL_OS":"ZL",
    "DYJetsMatchJ_OS":"ZJ",
    #"DYJetsSDB_OS":"Z_SDB",
    #"DYJetsQCD_OS":"ZQCD",
    #"DYJetsHMTSDB_SS":"Z_SS_HMT_SDB",
    "WJets_OS":"W",
    #"WJetsQCD_OS":"WQCD",
    #"WJetsQCDUp_OS":"WQCDUp",
    #"WJetsQCDDown_OS":"WQCDDown",
    #"WJetsQCDYield_OS":"WQCDYield",
    #"WJetsQCDYieldUp_OS":"WQCDYieldUp",
    #"WJetsQCDYieldDown_OS":"WQCDYieldDown",
    #"WJetsHigh_OS":"W_High",
    #"WJetsLow_OS":"W_Low",
    #"WJetsHMTSDB_SS":"W_SS_HMT_SDB",
    "TTbarMatchJ_OS":"TTJ",
    "TTbarMatchT_OS":"TTT",
    #"TTbar_OS":"TT",
    #"TTbarSDB_OS":"TT_SDB",
    #"TTbarHMTSDB_SS":"TT_SS_HMT_SDB",
    #"TTbarYield_OS":"TTYield",
    "ST_OS":"T",
    #"STQCD_OS":"TopQCD",
    "DiBosonMatchT_OS":"VVT",
    "DiBosonMatchJ_OS":"VVJ",
    #"DiBoson_OS":"VV",
    #"DiBosonQCD_OS":"VV_QCD",
    #"DiBosonSDB_OS":"VV_SDB",
    #"DiBosonHMTSDB_SS":"VV_SS_HMT_SDB",
    #"QCDEstimateHMTSDB_SS":"QCD_SS_HMT_SDB",
    "QCDEstimate":"QCD",
    "ggH120_OS":"ggH120",
    "qqH120_OS":"qqH120",
    "ggH125_OS":"ggH125",
    "qqH125_OS":"qqH125",
    "ggH130_OS":"ggH130",
    "qqH130_OS":"qqH130",
    "ZH120_OS":"ZH120",
    "ZH125_OS":"ZH125",
    "ZH130_OS":"ZH130",
    "WH120_OS":"WH120",
    "WH125_OS":"WH125",
    "WH130_OS":"WH130",
    "EWK2Jets_OS":"EWKZ",
    #"BkgErr":"BKGErr"
    }

#MT channel specific
#CMS_scale_t_13TeV and CMS_scale_j_13TeV are applied to all processes
#when another nuisance parameter is considered in the process, also histos Nuisance1_Nuisance2 should be added
nuisanceParams = {
    "CMS_htt_jetToTauFake_13TeV":(("",""),("ZJ","W","TTJ")),
    "CMS_htt_ZLShape_13TeV":(("mt_",""),("ZL","")),
    "CMS_htt_dyShape_13TeV":(("",""),("ZL","ZJ","ZTT")),
    "CMS_htt_ttbarShape_13TeV":(("",""),("TT","TTT","TTJ")),
    "QCDSFUncert_13TeV":(("mt_",""),("QCD","W")),
    "WSFUncert_13TeV":(("mt_",""),("W","")),
    "CMS_scale_gg_13TeV":(("",""),("ggH120","ggH125","ggH130"))
    }

nbins = {"mt_wjets_0jet_cr":(1,80,200),
         "mt_wjets_boosted_cr":(1,80,200),
         "mt_wjets_vbf_cr":(1,80,200),
         "mt_antiiso_0jet_cr":(4,40,200),
         "mt_antiiso_boosted_cr":(4,40,200),
         "mt_antiiso_vbf_cr":(4,40,200)
         }



categoryDirMade = True

for iCategory in xrange(0,len(categoryNames)):
    categoryName = categoryNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key
        if categoryName.count("antiiso")>0 and value!="QCD": 
            hName = hName+"noMuIso"
        hName = hName +"_"+str(iCategory)
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName, " is missing"
            histogram = TH1F(value,"",nbins[categoryName][0],nbins[categoryName][1],nbins[categoryName][2])
        histogram.SetName(value)
        histogram.Write()
        if value=="data_obs": continue

        for np in ("CMS_scale_j_13TeV", "CMS_scale_t_mt_13TeV"):
            histos = getSingleNPHistos(hName, np ,histogram)
            hUp = histos[0]
            hDown = histos[1]
            hUp.SetName(value+"_"+np+"Up")
            hUp.Write()
            hDown.SetName(value+"_"+np+"Down")
            hDown.Write()

        for key1, value1 in nuisanceParams.items():
            index = key1.find("13TeV")
            name = key1
            for chn in value1[0]:
                if (chn!=channel and chn!=""): continue
                name = key1[:index]+chn+key1[index:]
                if key1.count("Uncert")>0:
                  index = name.find("13TeV")
                  cat = categoryName.split("_")[2]
                  name = name[:index]+cat+"_"+name[index:]
                done1 = True

                for proc in value1[1]:
                  if proc!= value: continue
                  histos = getSingleNPHistos(hName, name, histogram)
                  hUp = histos[0]
                  hDown = histos[1]
                  hUp.SetName(value+"_"+name+"Up")
                  hUp.Write()
                  hDown.SetName(value+"_"+name+"Down")
                  hDown.Write()
                  '''
                  for np in ("CMS_scale_j_13TeV", "CMS_scale_t_mt_13TeV"):
                    histos = getDoubleNPHistos(hName, name, np ,histogram)
                    hUpUp = histos[0]
                    hUpDown = histos[1]
                    hDownUp = histos[2]
                    hDownDown = histos[3]
                    hUpUp.SetName(value+"_"+name+"Up_"+np+"Up")
                    hUpUp.Write()
                    hUpDown.SetName(value+"_"+name+"Up_"+np+"Down")
                    hUpDown.Write()
                    hDownUp.SetName(value+"_"+name+"Down_"+np+"Up")
                    hDownUp.Write()
                    hDownDown.SetName(value+"_"+name+"Down_"+np+"Down")
                    hDownDown.Write()
                  '''
                  if proc=="": break
                if (chn=="" or done1): break
