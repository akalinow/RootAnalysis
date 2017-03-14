from ROOT import TFile,TTree

input  = 'RootAnalysis_SynchNTupleTT_VBFMVAMET.root'
output = 'WAW_VBF_syncTree_tt_SM.root'
#input  = 'RootAnalysis_SynchNTupleMT_VBFMVAMET.root'
#output = 'WAW_VBF_syncTree_mt_SM.root'

fin = TFile(input)
t = fin.Get("Summary/tree")

fout = TFile(output,"recreate")
fout.mkdir("Summary")
fout.cd("Summary")
t.CloneTree().Write()

fout.Close()
fin.Close()
