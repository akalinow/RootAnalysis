from ROOT import TFile,TTree

input  = 'RootAnalysis_SynchMuTau.root'
output = 'WAW_SYNCFILE_SUSYGluGluToBBHToTauTau_M-1000_2016_mt_0.root'
#input  = 'RootAnalysis_SynchTT.root'
#output = 'WAW_SYNCFILE_SUSYGluGluToBBHToTauTau_M-1000_2016_tt_0.root'

fin = TFile(input)
t = fin.Get("Summary/tree")
#t = fin.Get("tree")

fout = TFile(output,"recreate")
fout.mkdir("Synch")
fout.cd("Synch")
t.CloneTree().Write()

fout.Close()
fin.Close()
