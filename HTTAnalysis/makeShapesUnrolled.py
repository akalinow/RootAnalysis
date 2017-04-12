#!/usr/bin/env python

from ROOT import *
import array
import numpy

WAW_fileName = "/cms/cms/akalinow/CMS/HiggsCP/Data/NTUPLES_28_03_2017/MT/Histograms/RootAnalysis_AnalysisMuTau.root"

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

categoryRootAnalysisNames = [
        "0jet", "boosted", "vbf",
        "0jet_W", "boosted_W", "vbf_W",
        "antiIso_0jet", "antiIso_boosted", "antiIso_vbf",
        ]

categoryCombineNames = list()

for i in xrange(0,len(categoryRootAnalysisNames)):
        tmp = categoryRootAnalysisNames[i]
        if tmp.count("_W"):
          tmp = tmp[:-2]
          tmp = "wjets_"+tmp
        tmp = channel+"_"+tmp
        if tmp.count("wjets")>0 or tmp.count("antiIso")>0:
                tmp = tmp + "_cr"
        tmp = tmp.replace('antiIso','antiiso')
        categoryCombineNames.append(tmp)

histogramsMap = {
    "Data":"data_obs",
    "DYJetsMatchT":"ZTT",
    "DYJetsMatchL":"ZL",
    "DYJetsMatchJ":"ZJ",
    "WJets":"W",
    "TTbarMatchJ":"TTJ",
    "TTbarMatchT":"TTT",
    "ST":"T",
    #"DiBosonMatchT":"VVT",
    #"DiBosonMatchJ":"VVJ",
    "DiBoson":"VV",
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

def rebinHisto(histo, categoryName):
    nbins = 0
    if categoryName.count("antiiso")>0:
        xbins = numpy.array([40.0,80.0,120.0,160.0,200.0])
        nbins = 4
    if categoryName.count("wjets")>0:
        xbins = numpy.array([80.0,200.0])
        nbins = 1
    if categoryName.count("qcd")>0:
        xbins = numpy.array([0.0,300.0])
        nbins = 1
    newHisto = histo.Rebin(nbins, "rebinned", xbins)
    return newHisto
    
def getHistogram(name, histo):
    if name.count("DiBoson")>0 and WAW_fileName.count("MuTau")>0: name=name.replace("DiBoson","DiBosonMatchT")
    h = WAW_file.Get(name)
    txtFile=open("histogramSearch.txt","a")
    if h==None :
        print "Missing histogram: ",name
        txtFile.write("Search for histo: "+name+": 0\n")
        h = histo.Clone()
    else: 
        txtFile.write("Search for histo: "+name+": 1\n")
        
    if name.count("DiBoson")>0 and WAW_fileName.count("MuTau")>0:
        name=name.replace("MatchT","MatchJ")
        hJ = WAW_file.Get(name)
        if hJ==None :
          print "Missing histogram: ",name
          txtFile.write("Search for histo: "+name+": 0\n")
          hJ = histo.Clone()
        else:
          txtFile.write("Search for histo: "+name+": 1\n")
          h.Add(hJ)
    txtFile.close()
    return h

def getSingleNPHistos(prefix, np, histo):
        
    hName = prefix+"_"+np
    hUp = getHistogram(hName+"Up",histo)
    hDown = getHistogram(hName+"Down",histo)
    return (hUp, hDown)

#basic categories
categoryDirMade = False

hData = 0

open("histogramSearchMT.txt","w").write("Looking for histograms for Combine \n")

print "MAIN REGION\n\n\n\n\n"

for iCategory in xrange(0,len(categoryCombineNames)):
    categoryName = categoryCombineNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key+"_"+categoryRootAnalysisNames[iCategory]
        templateHisto = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram = getHistogram(hName, templateHisto)
        histogram.SetName(value)
        histogram.Write()

        if value=="data_obs": continue

        for nuisanceParam in nuisanceParams:
            nuisanceParam = nuisanceParam.replace("CHANNEL",channel)
            cat = categoryName
            if cat.count("0jet")>0: cat = "0jet"
            elif cat.count("boosted")>0: cat = "boosted"
            elif cat.count("vbf")>0: cat = "vbf"
            nuisanceParam = nuisanceParam.replace("CAT",cat)
            
            histos = getSingleNPHistos(histoPrefix[categoryName] + key+"_"+categoryRootAnalysisNames[iCategory], nuisanceParam, histogram)
            if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0:  nuisanceParam=nuisanceParam.replace("vbf","VBF")
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
    #"DiBosonMatchT":"VVT",
    #"DiBosonMatchJ":"VVJ",
    "DiBoson":"VV",
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
    "CMS_htt_ZLShape_CHANNEL_13TeV":("ZL",),
    "CMS_htt_dyShape_13TeV":("ZL","ZJ","ZTT"),
    "CMS_htt_ttbarShape_13TeV":("TT","TTT","TTJ"),
    "QCDSFUncert_CHANNEL_CAT_13TeV":("QCD","W"),
    "WSFUncert_CHANNEL_CAT_13TeV":("W",),
    "CMS_scale_gg_13TeV":("ggH120","ggH125","ggH130"),
    "CMS_htt_zmumuShape_CAT_13TeV":("ZTT","ZL", "ZJ")
    }
    
nbins = {"mt_wjets_0jet_cr":(1,80,200),
         "mt_wjets_boosted_cr":(1,80,200),
         "mt_wjets_vbf_cr":(1,80,200),
         "mt_antiiso_0jet_cr":(4,40,200),
         "mt_antiiso_boosted_cr":(4,40,200),
         "mt_antiiso_vbf_cr":(4,40,200)
         }



categoryDirMade = True

for iCategory in xrange(0,len(categoryCombineNames)):
    categoryName = categoryCombineNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key
        hName = hName +"_"+categoryRootAnalysisNames[iCategory]
        templateHisto = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram = getHistogram(hName, templateHisto)
        histogram=rebinHisto(histogram, categoryName)
        histogram.SetName(value)
            
        histogram.Write()
        if value=="data_obs": continue

        for np in ("CMS_scale_j_13TeV", "CMS_scale_t_mt_13TeV"):
            histos = getSingleNPHistos(hName, np ,histogram)
            hUp = histos[0]
            hDown = histos[1]
            hUp=rebinHisto(hUp, categoryName)
            hDown=rebinHisto(hDown, categoryName)
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
            nuisanceParamRA = nuisanceParam.replace("CAT",categoryRootAnalysisNames[iCategory])
            nuisanceParam = nuisanceParam.replace("CAT",cat)

            for proc in value1:
                  if proc!= value: continue
                  histos = getSingleNPHistos(hName, nuisanceParamRA, histogram)
                  if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0:
                    nuisanceParam = nuisanceParam.replace("vbf","VBF")
                  hUp = histos[0]
                  hDown = histos[1]
                  hUp=rebinHisto(hUp, categoryName)
                  hDown=rebinHisto(hDown, categoryName)
                  hUp.SetName(value+"_"+nuisanceParam+"Up")
                  hUp.Write()
                  hDown.SetName(value+"_"+nuisanceParam+"Down")
                  hDown.Write()
                

##########################################################
######### tauH-tauH channel###############################
##########################################################


WAW_fileName = "RootAnalysis_AnalysisTT.root"

channel="tt"

SYNCH_fileName = "htt_tt.inputs-sm-13TeV-2D.root"

WAW_file = TFile(WAW_fileName)
SYNCH_file = TFile(SYNCH_fileName,"RECREATE")

#WAW to SYNCH histograms names map
histoPrefix = {
         "tt_0jet":"HTTAnalyzer/h1DMassSV",
         "tt_boosted":"HTTAnalyzer/h1DUnRollHiggsPtMassSV",
         "tt_vbf":"HTTAnalyzer/h1DUnRollMjjMassSV"
         }

#define number of bins in each category (nbinsX, nbinsY) in case when you need to create an empty histo
nbins = {
         "tt_0jet":(30,1),
         "tt_boosted":(12,4),
         "tt_vbf":(12,4)
         }

categoryRootAnalysisNames = [
        "0jet", "boosted", "vbf",
        "0jet_QCD", "boosted_QCD", "vbf_QCD"
        ]

categoryCombineNames = list()

for i in xrange(0,len(categoryRootAnalysisNames)):
        tmp = categoryRootAnalysisNames[i]
        if tmp.count("_QCD"):
          tmp = tmp.replace("QCD","qcd")
          tmp+="_cr"
        tmp = channel+"_"+tmp
        categoryCombineNames.append(tmp)

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
    #"DiBoson":"VV",
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

#basic categories
categoryDirMade = False

hData = 0

open("histogramSearchMT.txt","w").write("Looking for histograms for Combine \n")

print "MAIN REGION\n\n\n\n\n"

for iCategory in xrange(0,len(categoryCombineNames)):
    categoryName = categoryCombineNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key+"_"+categoryRootAnalysisNames[iCategory]
        templateHisto = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram = getHistogram(hName, templateHisto)
        histogram.SetName(value)
        histogram.Write()

        if value=="data_obs": continue

        for nuisanceParam in nuisanceParams:
            nuisanceParam = nuisanceParam.replace("CHANNEL","mt")
            cat = categoryName
            if cat.count("0jet")>0: cat = "0jet"
            elif cat.count("boosted")>0: cat = "boosted"
            elif cat.count("vbf")>0: cat = "vbf"
            nuisanceParam = nuisanceParam.replace("CAT",cat)
            
            histos = getSingleNPHistos(histoPrefix[categoryName] + key+"_"+categoryRootAnalysisNames[iCategory], nuisanceParam, histogram)
            if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0:  nuisanceParam=nuisanceParam.replace("vbf","VBF")
            nuisanceParam=nuisanceParam.replace("mt","tt")
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
         "tt_0jet_qcd_cr":"HTTAnalyzer/h1DMassVis",
         "tt_boosted_qcd_cr":"HTTAnalyzer/h1DMassVis",
         "tt_vbf_qcd_cr":"HTTAnalyzer/h1DMassVis"
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

#TT channel specific
#CMS_scale_t_13TeV and CMS_scale_j_13TeV are applied to all processes
#when another nuisance parameter is considered in the process, also histos Nuisance1_Nuisance2 should be added
nuisanceParams = {
    "CMS_htt_jetToTauFake_13TeV":("ZJ","W","TTJ", "VVJ"),
    "CMS_htt_ZLShape_CHANNEL_13TeV":("ZL",),
    "CMS_htt_dyShape_13TeV":("ZL","ZJ","ZTT"),
    "CMS_htt_ttbarShape_13TeV":("TT","TTT","TTJ"),
    "CMS_scale_gg_13TeV":("ggH120","ggH125","ggH130"),
    "CMS_htt_zmumuShape_CAT_13TeV":("ZTT","ZL", "ZJ", "EWKZ")
    }
    
nbins = {"tt_0jet_qcd_cr":(1,1),
         "tt_boosted_qcd_cr":(1,1),
         "tt_vbf_qcd_cr":(1,1)
         }



categoryDirMade = True

for iCategory in xrange(0,len(categoryCombineNames)):
    categoryName = categoryCombineNames[iCategory]
    if categoryName not in histoPrefix.keys(): continue
    if categoryDirMade: gDirectory.cd("..")
    gDirectory.mkdir(categoryName)
    gDirectory.cd(categoryName)
    categoryDirMade=True

    for key,value in histogramsMap.iteritems():
        hName = histoPrefix[categoryName] + key
        hName = hName +"_"+categoryRootAnalysisNames[iCategory]
        templateHisto = TH1F(value,"",nbins[categoryName][0]*nbins[categoryName][1],0.5,nbins[categoryName][0]*nbins[categoryName][1] + 0.5)
        histogram = getHistogram(hName, templateHisto)
        histogram=rebinHisto(histogram, categoryName)
        histogram.SetName(value)
            
        histogram.Write()
        if value=="data_obs": continue

        for np in ("CMS_scale_j_13TeV", "CMS_scale_t_mt_13TeV"):
            histos = getSingleNPHistos(hName, np ,histogram)
            hUp = histos[0]
            hDown = histos[1]
            hUp=rebinHisto(hUp, categoryName)
            hDown=rebinHisto(hDown, categoryName)
            np=np.replace("mt","tt")
            hUp.SetName(value+"_"+np+"Up")
            hUp.Write()
            hDown.SetName(value+"_"+np+"Down")
            hDown.Write()

        for nuisanceParam, value1 in nuisanceParams.items():
            nuisanceParam = nuisanceParam.replace("CHANNEL","mt")
            cat = categoryName
            if cat.count("0jet")>0: cat = "0jet"
            elif cat.count("boosted")>0: cat = "boosted"
            elif cat.count("vbf")>0: cat = "vbf"
            nuisanceParamRA = nuisanceParam.replace("CAT",categoryRootAnalysisNames[iCategory])
            nuisanceParam = nuisanceParam.replace("CAT",cat)
            nuisanceParam = nuisanceParam.replace("mt","tt")

            for proc in value1:
                  if proc!= value: continue
                  histos = getSingleNPHistos(hName, nuisanceParamRA, histogram)
                  if nuisanceParam.count("zmumuShape")>0 and cat.count("vbf")>0:
                    nuisanceParam = nuisanceParam.replace("vbf","VBF")
                  hUp = histos[0]
                  hDown = histos[1]
                  hUp=rebinHisto(hUp, categoryName)
                  hDown=rebinHisto(hDown, categoryName)
                  hUp.SetName(value+"_"+nuisanceParam+"Up")
                  hUp.Write()
                  hDown.SetName(value+"_"+nuisanceParam+"Down")
                  hDown.Write()
                

