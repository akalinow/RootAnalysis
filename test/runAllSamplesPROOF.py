#!/usr/bin/env python

import os,re,time
import getSelectedEventList


from VBFZtautau_mutau_cfg import *
#########################################
def processSingleSample(paths, sampleName, sampleWeight, HLTprocessName,isMC=True,useFilter=False):
	process.source.fileNames = cms.untracked.vstring()
	for path in paths:
		for fname in os.listdir(path):
			if re.search(".root",fname)!=None and  os.path.getmtime(path+fname)-time.time()<-60:
				process.source.fileNames.append('file:'+path+fname)
	process.configurationMetadata.name=sampleName
	process.configurationMetadata.preselectionEff=sampleWeight
	process.TriggerAnalyzer.triggerResultsLabel = cms.InputTag("TriggerResults","",HLTprocessName)
	if isMC:
		process.DiTau.useEventWeight = cms.untracked.bool(True)
		process.DiTau.useMCTruthData = cms.bool(True)
		process.DiTau.matchTriggerFilterNames = cms.vstring("hltSingleMuIsoL3IsoFiltered15",
								    "hltOverlapFilterIsoMu12IsoPFTau10")
		process.TriggerAnalyzer.triggerItemNames = cms.vstring("HLT_IsoMu15_v5","HLT_IsoMu12_LooseIsoPFTau10_v2")
	else:
		process.DiTau.useEventWeight = cms.untracked.bool(False)
		process.DiTau.useMCTruthData = cms.bool(False)
		process.TriggerAnalyzer.triggerItemNames = cms.vstring("HLT_IsoMu17_v5",
								       "HLT_IsoMu17_v6",
								       "HLT_IsoMu17_v7",
								       "HLT_IsoMu17_v8",
								       "HLT_IsoMu17_v9",
								       "HLT_IsoMu17_v10",
								       "HLT_IsoMu17_v11",
								       "HLT_IsoPFTau35_Trk20_MET45_v1",
								       "HLT_IsoPFTau35_Trk20_MET45_v2",
								       "HLT_IsoPFTau35_Trk20_MET45_v4",
                                                                       "HLT_IsoMu15_LooseIsoPFTau15_v2",
                                                                       "HLT_IsoMu15_LooseIsoPFTau15_v4",
                                                                       "HLT_IsoMu15_LooseIsoPFTau15_v6",
			                                               "HLT_IsoMu15_LooseIsoPFTau20_v1",
								       "HLT_IsoMu15_LooseIsoPFTau20_v2",
								       "HLT_IsoMu15_LooseIsoPFTau20_v4",
								       "HLT_IsoMu12_LooseIsoPFTau10_v1",
								       "HLT_IsoMu12_LooseIsoPFTau10_v2",
								       "HLT_IsoMu12_LooseIsoPFTau10_v4")
		process.DiTau.matchTriggerFilterNames = cms.vstring("hltSingleMuIsoL3IsoFiltered17",
								    "hltOverlapFilterIsoMu12IsoPFTau10",
								    "hltOverlapFilterIsoMu15IsoPFTau15",
								    "hltOverlapFilterIsoMu15IsoPFTau20")
	###
	path = "/scratch/scratch0/akalinow/3500GeV/CMSSW_4_2_4_patch1/src/PFAnalyses/VBFHTauTau/test/06_07_2011/"
	fileName = path+"PFAnalysis_WJetsToLNu.root"
	flavour = "HPSSingleMu"
	if useFilter:
		getSelectedEventList.getEventsToProcess(fileName,flavour)
	###	
	out = open('tmpConfig.py','w')
	out.write(process.dumpPython())
	out.close()
	os.system("taskset -c 1,2,3,4,5,6,7 root -q -b runProof.C")
	#os.system("rm -f tmpConfig.py")
#########################################
#The weight is gives in following components:
# BR(2tau->mu + jet)*sigma*BR(X->tau tau)*(eff gen preselection)*(eff reco preselection)*(conversion pb->fb)
#
# Note1: In the case of Z->tau tau the BR(X->tau tau) is included in the cross section
# Note2: In the case of W->tau the BR(X->tau) is included in the cross section
# Note3: In the case of tt thre is no decay forcing
# Note4: eff gen preselection: is efficiency to have a generator level muon/tau in the acceptance area
# Note5: eff reco preselection: is efficiency to have a reco level DiTau
# Note6: Cross section taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
#        With lepton universality assumed.
#########################################
rootPath = "/scratch/scratch0/akalinow/3500GeV/CMSSW_4_2_4_patch1/src/PFAnalyses/VBFHTauTau/test/Crab/"

################
dataSingleMuReReco = {
	"/SingleMu/Run2011A-May10ReReco-v1/AOD":\
	("Run2011A_May10ReReco-v1_SingleMu",-1.0,"HLT",rootPath+"SingleMu/Run2011A-May10ReReco-v1/AOD/crab_0_110707_013729/res/")
	}
################
dataSingleMuPromptReco = {
	"/SingleMu/Run2011A-PromptReco-v4/AOD":\
	("Run2011A-PromptReco-v4_SingleMu",-1.0,"HLT",rootPath+"SingleMu/Run2011A-PromptReco-v4/AOD/crab_0_110707_013541/res/")
	}
################
DY = 	{"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Spring11-PU_S1_START311_V1G1-v1/AODSIM/":\
	 ("DYJetsToLL",-3048.0,"HLT",rootPath+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM/crab_0_110701_095642/res/"),
	 }
################
WJets = {"WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM/":\
	("WJetsToLNu",-31314,"HLT",rootPath+"/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM/crab_0_110701_095802/res/"),
	 }
################
datasets = dataSingleMuReReco
#datasets.update(dataSingleMuPromptReco)

#datasets = dataSingleMuPromptReco

for key in datasets.keys():
		aDataSet = key
		sampleName = datasets[key][0]
		scaleFactor = datasets[key][1]
		HLTProcessName = datasets[key][2]
		samplePath = datasets[key][3]
		processSingleSample([samplePath],sampleName,scaleFactor,HLTProcessName,False,False)

'''
datasets = DY
datasets.update(WJets)

for key in datasets.keys():
		aDataSet = key
		sampleName = datasets[key][0]
		scaleFactor = datasets[key][1]
		HLTProcessName = datasets[key][2]
		samplePath = datasets[key][3]
		processSingleSample([samplePath],sampleName,scaleFactor,HLTProcessName,True,False)
'''
