import read_tree as rt
import pickle
import sys

def root2pickle(root_file = "../build/RootAnalysis_AnalysisMuTau.root",
								tree_path="Summary/tree",
								out_file="htt_features.pkl"):
	"""Imports data from the .root file using read_tree and saves them to pickle file
	data can be then accessed with:
	```
	with open("htt_features.pkl", "rb") as input:
		legs, jets, global_params, properties = pickle.load(input)
	```
	"""
	legs, jets, global_params, properties = rt.read_tree(root_file, tree_path)
	with open("htt_features.pkl", "wb") as output:
			pickle.dump((legs,jets, global_params, properties), output)
	print("Data from " + root_file + " saved to " + out_file)

if __name__ == "__main__":
	root2pickle(*sys.argv[1:])
		
	
