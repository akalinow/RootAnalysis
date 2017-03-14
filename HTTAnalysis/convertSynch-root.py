from ROOT import TFile,TTree

#input  = 'RootAnalysis_SynchMuTau.root'
#output = 'WAW_SYNCFILE_SUSYGluGluToBBHToTauTau_M-1000_2016_mt_0.root'
input  = 'RootAnalysis_SynchTT.root'
output = 'WAW_SYNCFILE_SUSYGluGluToBBHToTauTau_M-1000_2016_tt_0.root'

fin = TFile(input)
t = fin.Get("Summary/tree")

fout = TFile(output,"recreate")
fout.mkdir("Summary")
fout.cd("Summary")
t.CloneTree().Write()

fout.Close()
fin.Close()
