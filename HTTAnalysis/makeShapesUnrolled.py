#!/usr/bin/env python

from ROOT import *
import array
import numpy

WAW_fileName = "RootAnalysis_AnalysisMuTau_Artur.root"

channel="mt"

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

categoryEnum = [
        "jet0_low", "jet0_high",
        "jet1_low", "jet1_high",
        "vbf_low", "vbf_high",
        "jet0", "boosted", "vbf",
        "wjets_jet0", "wjets_boosted", "wjets_vbf",
        "wjets_qcd_jet0", "wjets_qcd_boosted", "wjets_qcd_vbf",
        "qcd_jet0", "qcd_boosted", "qcd_vbf",
        "qcd_ss_jet0", "qcd_ss_boosted", "qcd_ss_vbf",
        "ss_jet0", "ss_boosted", "ss_vbf",
        "antiIso_jet0", "antiIso_boosted", "antiIso_vbf",
        "antiIso_qcd_jet0", "antiIso_qcd_boosted", "antiIso_qcd_vbf"
        ]

categoryNames = list()

for i in xrange(0,len(categoryEnum)):
        tmp = categoryEnum[i]
        tmp = channel+"_"+tmp
        if tmp.count("wjets")>0 or tmp.count("antiIso")>0:
                tmp = tmp + "_cr"
        tmp = tmp.replace('jet0','0jet')
        tmp = tmp.replace('antiIso','antiiso')
        categoryNames.append(tmp)

histogramsMap = {
    "Data":"data_obs",
    "DYJetsMatchT":"ZTT",
    "DYJetsMatchL":"ZL",
    "DYJetsMatchJ":"ZJ",
    "WJets":"W",
    "TTbarMatchJ":"TTJ",
    "TTbarMatchT":"TTT",
    "ST":"T",
    "DiBosonMatchT":"VVT",
    "DiBosonMatchJ":"VVJ",
    "QCDEstimate":"QCD",
    "ggHTT120":"ggH120",
    "qqHTT120":"qqH120",
    "ggHTT125":"ggH125",
    "qqHTT125":"qqH125",
    "ggHTT130":"ggH130",
    "qqHTT130":"qqH130",
    "ZHTT120":"ZH120",
    "WHTT120":"WH120",
    "ZHTT125":"ZH125",
    "WHTT125":"WH125",
    "ZHTT130":"ZH130",
    "WHTT130":"WH130",
    "EWK2Jets":"EWKZ"
    }

#according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#Systematic_uncertainties
nuisanceParams = [
    #"CMS_shape_t_CHANNEL_13TeV",
    "CMS_scale_t_CHANNEL_13TeV",
    #"CMS_scale_e_CHANNEL_13TeV",
    "CMS_scale_j_13TeV",
    "CMS_htt_jetToTauFake_13TeV",
    "CMS_htt_ZLShape_CHANNEL_13TeV",
    "CMS_htt_dyShape_13TeV",
    "CMS_htt_ttbarShape_13TeV",
    "QCDSFUncert_CHANNEL_CAT_13TeV",
    "WSFUncert_CHANNEL_CAT_13TeV",
    "CMS_scale_gg_13TeV",
    "CMS_htt_zmumuShape_CAT_13TeV",
    ]
#AAAAAAAAAAAAAAa brakuje CMS_ przy scale_gg i nizej tez
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

def rebinHisto(histo, categoryName):
    if categoryName.count("antiiso")>0:
        histo.Rebin(4)
        newHisto = TH1F(histo.GetName(),histo.GetTitle(),4,40,200)
        for i in xrange(1,5):
            newHisto.SetBinContent(i, histo.GetBinContent(i+1))
        return newHisto
    if categoryName.count("wjets")>0:
        newHisto = TH1F(histo.GetName(),histo.GetTitle(),1,80,200)
        newHisto.SetBinContent(1, histo.Integral(histo.GetBin(81),histo.GetBin(199)))
        return newHisto

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
        hName = histoPrefix[categoryName] + key
        histogram = WAW_file.Get(hName)
        if(histogram==None):
            print hName,"is missing"
            histogram = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram.SetName(value)
        histogram.Write()

        if value=="data_obs": continue

        for nuisanceParam in nuisanceParams:
            nuisanceParam = nuisanceParam.replace("CHANNEL",channel)
            cat = categoryName
            if cat.count("0jet")>0: cat = "0jet"
            elif cat.count("boosted")>0: cat = "boosted"
            elif cat.count("vbf")>0: cat = "vbf"
            if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0: cat = "VBF"
            nuisanceParam = nuisanceParam.replace("CAT",cat)
            
            histos = getSingleNPHistos(histoPrefix[categoryName] + key, nuisanceParam, histogram)
            histogramUp = histos[0]
            histogramUp.SetName(value+"_"+nuisanceParam+"Up")
            histogramUp.Write()

            histogramDown = histos[1]
            histogramDown.SetName(value+"_"+nuisanceParam+"Down")
            histogramDown.Write()

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
    "Data":"data_obs",
    "DYJetsMatchT":"ZTT",
    "DYJetsMatchL":"ZL",
    "DYJetsMatchJ":"ZJ",
    #"DYJetsSDB":"Z_SDB",
    #"DYJetsQCD":"ZQCD",
    #"DYJetsHMTSDB_SS":"Z_SS_HMT_SDB",
    "WJets":"W",
    #"WJetsQCD":"WQCD",
    #"WJetsQCDUp":"WQCDUp",
    #"WJetsQCDDown":"WQCDDown",
    #"WJetsQCDYield":"WQCDYield",
    #"WJetsQCDYieldUp":"WQCDYieldUp",
    #"WJetsQCDYieldDown":"WQCDYieldDown",
    #"WJetsHigh":"W_High",
    #"WJetsLow":"W_Low",
    #"WJetsHMTSDB_SS":"W_SS_HMT_SDB",
    "TTbarMatchJ":"TTJ",
    "TTbarMatchT":"TTT",
    #"TTbar":"TT",
    #"TTbarSDB":"TT_SDB",
    #"TTbarHMTSDB_SS":"TT_SS_HMT_SDB",
    #"TTbarYield":"TTYield",
    "ST":"T",
    #"STQCD":"TopQCD",
    "DiBosonMatchT":"VVT",
    "DiBosonMatchJ":"VVJ",
    #"DiBoson":"VV",
    #"DiBosonQCD":"VV_QCD",
    #"DiBosonSDB":"VV_SDB",
    #"DiBosonHMTSDB_SS":"VV_SS_HMT_SDB",
    #"QCDEstimateHMTSDB_SS":"QCD_SS_HMT_SDB",
    "QCDEstimate":"QCD",
    "ggHTT120":"ggH120",
    "qqHTT120":"qqH120",
    "ggHTT125":"ggH125",
    "qqHTT125":"qqH125",
    "ggHTT130":"ggH130",
    "qqHTT130":"qqH130",
    "ZHTT120":"ZH120",
    "ZHTT125":"ZH125",
    "ZHTT130":"ZH130",
    "WHTT120":"WH120",
    "WHTT125":"WH125",
    "WHTT130":"WH130",
    "EWK2Jets":"EWKZ",
    #"BkgErr":"BKGErr"
    }

#MT channel specific
#CMS_scale_t_13TeV and CMS_scale_j_13TeV are applied to all processes
#when another nuisance parameter is considered in the process, also histos Nuisance1_Nuisance2 should be added
nuisanceParams = {
    "CMS_htt_jetToTauFake_13TeV":("ZJ","W","TTJ"),
    "CMS_htt_ZLShape_CHANNEL_13TeV":("ZL",""),
    "CMS_htt_dyShape_13TeV":("ZL","ZJ","ZTT"),
    "CMS_htt_ttbarShape_13TeV":("TT","TTT","TTJ"),
    "QCDSFUncert_CHANNEL_CAT_13TeV":("QCD","W"),
    "WSFUncert_CHANNEL_CAT_13TeV":("W",""),
    "CMS_scale_gg_13TeV":("ggH120","ggH125","ggH130")
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

        for nuisanceParam, value1 in nuisanceParams.items():
            nuisanceParam = nuisanceParam.replace("CHANNEL",channel)
            cat = categoryName
            if cat.count("0jet")>0: cat = "0jet"
            elif cat.count("boosted")>0: cat = "boosted"
            elif cat.count("vbf")>0: cat = "vbf"
            if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0: cat = "VBF"
            nuisanceParam = nuisanceParam.replace("CAT",cat)

            for proc in value1:
                  if proc!= value: continue
                  histos = getSingleNPHistos(hName, nuisanceParam, histogram)
                  hUp = histos[0]
                  hDown = histos[1]
                  hUp.SetName(value+"_"+nuisanceParam+"Up")
                  hUp.Write()
                  hDown.SetName(value+"_"+nuisanceParam+"Down")
                  hDown.Write()
                
