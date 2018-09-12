#!/usr/bin/env python
# -*- coding: utf-8 -*-

# \brief Script to read data from TTree inside a .root file into python lists.
# \author Rafal Maselek

from ROOT import TTree, TFile, TObjArray
import matplotlib.pyplot as plt

def read_tree(file_name = "../data/dummy.root", tree_name = "Summary/tree"):
	"""
	Function to read from a TTree inside a root file. Compatible with MLAnalyzer output.
	:param file_name: Name/path of the input .root file with the TTree inside.
	:param tree_name: Name/path to the TTree inside the .root file.
	:return: legs, jets, global_params, properties
		legs -- list of lists representing p4 of all legs in analysis; each list contains 4 sub-lists representing E,
		pX, pY, pZ; each value corresponds to single event
		jets -- list of lists representing p4 of all jets in analysis; each list contains 4 sub-lists representing E,
		pX, pY, pZ; each value corresponds to single event
		global_params -- dictionary; each key corresponds to a global parameter, each val is a list containing values
		of concrete global parameter; each value corresponds to single event
		properties -- dictionary; each key corresponds to a particle property (it includes particle identifier),
		each val is a list containing values of concrete property; each value corresponds to single event
	"""

	# Disk work
	print("[ML]\tOpening root file for python conversion.")
	file = TFile.Open(file_name, 'read')
	tree = file.Get(tree_name)
	leaves = tree.GetListOfLeaves()

	# define dynamically a python class containing root Leaves objects
	class PyListOfLeaves(dict) :
	    pass

	# create an istance
	pyl = PyListOfLeaves()

	# Deal with leaves (branches) by putting them to pyl
	leaves_names = []
	for i in range(0, leaves.GetEntriesFast()):
	    leaf = leaves.At(i)
	    name = leaf.GetName()
	    leaves_names.append(name)
	    # add dynamically attribute to my class 
	    pyl.__setattr__(name, leaf)

	# containers for loaded data
	jets = []
	legs = []
	global_params_names = []
	properties_names = []

	# endings of leaves that correspond to four-momentum
	p4 = ["E", "pX", "pY", "pZ"]
	leaves_names.sort() # IMPORTANT! Names are now sorted!

	# Code below is responsible for preparing containers to match the data provided
	last_no = '0'
	for name in leaves_names:
		sub = name.split('_')
		# print(sub)
		if len(sub) >= 3:
			if sub[2] in p4:
				# We're dealing with four-momentum
				if sub[0] == 'leg':
					if sub[1] != last_no:
						last_no = sub[1]
						legs.append([]) # append list for new leg
						legs[-1].append([]) # append list for E
					else:
						legs[-1].append([]) # append list for pX,pY,pZ
				elif sub[0] == 'jet':
					if sub[1] != last_no:
						last_no = sub[1]
						jets.append([]) # append list for new jet
						jets[-1].append([]) # append list for E
					else:
						jets[-1].append([]) # append list for pX,pY,pZ
				else:
					raise Exception("[ERROR] P4 FOR UNKNOWN PARTICLE!")
			else:
				# We're dealing with particle parameters
				properties_names.append(name) # we use full name because may be for jets and legs

		elif sub[0] != 'eventWeight':
			# We're dealing with global parameters. We discard eventWeight
			global_params_names.append(sub[0])

	# Assign dictionaries for properties and global parameters
	global_params = {key: [] for key in global_params_names}
	properties = {key: [] for key in properties_names}
		

	# Loop over all entries and put data into proper containers
	print("[ML]\tReading data from TTree.")
	nev = tree.GetEntries()
	for iev in range(0, nev) :
	    tree.GetEntry(iev)
	    # get values from the tree using Python class pyl which contains leaves objects
	    for ll in range(0, len(legs)):
	    	for ii, entry in enumerate(p4):
	    		legs[ll][ii].append(pyl.__getattribute__("leg_"+str(ll+1)+"_"+entry).GetValue()) # append p4 values
	    for jj in range(0, len(jets)):
	    	for ii, entry in enumerate(p4):
	    		jets[jj][ii].append(pyl.__getattribute__("jet_"+str(jj+1)+"_"+entry).GetValue()) # append p4 values
	    for g_param in global_params_names:
	    	global_params[g_param].append(pyl.__getattribute__(g_param).GetValue())
	    for prop in properties_names:
	    	properties[prop].append(pyl.__getattribute__(prop).GetValue())

	###### optional check (uncomment what is below)
	# print(properties)
	# plt.hist(properties[properties_names[0] ])
	# plt.show()

	print("[ML]\tConversion to python successful!")
	return legs, jets, global_params, properties

if __name__ == "__main__":
    read_tree()
