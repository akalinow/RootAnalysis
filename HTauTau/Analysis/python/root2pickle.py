import read_tree as rt
import pickle
import sys
import argparse
import os

def root2pickle(root_file = "./RootAnalysis_SVfitAnalysisMuTau.root",
                out_file = "htt_features.pkl",
	        tree_path = "Summary/tree"):
								
	"""Imports data from the .root file using read_tree and saves them to pickle file
	data can be then accessed with:
	```
	with open("htt_features.pkl", "rb") as input:
		legs, jets, global_params, properties = pickle.load(input)
	```
	"""
	legs, jets, global_params, properties = rt.read_tree(root_file, tree_path)
	with open(out_file, "wb") as output:
			pickle.dump((legs,jets, global_params, properties), output)
	print("Data from " + root_file + " saved to " + out_file)

if __name__ == "__main__":

        parser = argparse.ArgumentParser()

        parser.add_argument('--input', type=str,
                            default=os.path.join(os.getenv('PWD', './'),
                                                 'RootAnalysis_SVfitAnalysisMuTau.root'),
                            help='Input ROOT file.')

        parser.add_argument('--output', type=str,
                            default=os.path.join(os.getenv('PWD', './'),
                                                 'htt_features.pkl'),
                            help='Output Python pickled file.')

        parser.add_argument('--tree', type=str,
                            default='Summary/tree',
                            help='Full path to the TTree inside the input file.')

        FLAGS, unparsed = parser.parse_known_args()
	root2pickle(FLAGS.input, FLAGS.output, FLAGS.tree)
		
	
