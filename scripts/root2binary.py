#!/usr/bin/env python
# -*- coding: utf-8 -*-

# \brief Example script that converse root data to binary tf.
# \author Rafal Maselek

import read_tree as rt
import io_tf_binary as io
import sys
import tensorflow as tf


# template for reader/writer
binary_template = io.Io_tf_binary('out_file',
			{'leg1_p4':(4,'f'), 
			'leg2_p4':(4,'f'),
			'jet1_p4':(4,'f'),
			'jet2_p4':(4,'f'),
			'nJets30':(1,'i'),
			'higgsPT':(1,'f'),
			'BJetBetaScore':(1,'f'),
			'visMass':(1,'f'),
			'higgsMassTrans':(1,'f'),
			'leg_2_charge':(1,'f'),
			'leg_2_byCombinedIsolationDeltaBetaCorrRaw3Hits':(1,'f'),
			'leg_1_combreliso':(1,'f'),
			'leg_2_byIsolationMVArun2v1DBoldDMwLTraw':(1,'f'),
			'leg_2_decayMode':(1,'f'),
			'leg_1_charge':(1,'f')
			})

# define  generator
def HTT_generator(legs, jets, int_pars, float_pars):
    for iev in range(len(legs[0][0])):
		# unpacking data
		leg1_p4 = [legs[0][ii][iev] for ii in range(4)]
		leg2_p4 = [legs[1][ii][iev] for ii in range(4)]
		jet1_p4 = [jets[0][ii][iev] for ii in range(4)]
		jet2_p4 = [jets[1][ii][iev] for ii in range(4)]
		
    	# TODO: Is it possible to be more generic?
		yield {'leg1_p4':leg1_p4, 
				'leg2_p4':leg2_p4,
        		'jet1_p4':jet1_p4,
        		'jet2_p4':jet2_p4,
        		'nJets30':int_pars['nJets30'][iev],
        		'higgsPT':float_pars['higgsPT'][iev],
        		'BJetBetaScore':float_pars['BJetBetaScore'][iev],
        		'visMass':float_pars['visMass'][iev],
        		'higgsMassTrans':float_pars['higgsMassTrans'][iev],
        		'leg_2_charge':float_pars['leg_2_charge'][iev],
        		'leg_2_byCombinedIsolationDeltaBetaCorrRaw3Hits':float_pars['leg_2_byCombinedIsolationDeltaBetaCorrRaw3Hits'][iev],
        		'leg_1_combreliso':float_pars['leg_1_combreliso'][iev],
        		'leg_2_byIsolationMVArun2v1DBoldDMwLTraw':float_pars['leg_2_byIsolationMVArun2v1DBoldDMwLTraw'][iev],
        		'leg_2_decayMode':float_pars['leg_2_decayMode'][iev],
        		'leg_1_charge':float_pars['leg_1_charge'][iev]
        		},0 #@Pawel Czy my nie mamy uczyc bez nadzoru? Wstawiam tu 0...


def save2binary(in_file="../data/dummy.root", tree_path="Summary/tree", out_file="output/example"):

	# load data from root file
	legs, jets, global_params, properties = rt.read_tree(in_file, tree_path)

	print("[ML]\tNumber of legs: {}".format(len(legs)))
	print("[ML]\tNumber of jets: {}".format(len(jets)))
	print("[ML]\tGlobal parameters from data:")
	for key, value in global_params.iteritems():
		print("\t*  "+key)
	print("[ML]\tParticle properties from data:")
	for key, value in properties.iteritems():
		print("\t*  "+key)

	# do some stuf with data
	# stuff
	# more stuff
	
	ints = {}
	floats = {}

	if "nJets30" in global_params:
		ints = {"nJets30": [int(x) for x in global_params["nJets30"]]} # conversion to integer

	for key in global_params:
		if key != "nJets30":
			floats[key] = global_params[key]
	for key in properties:
		floats[key] = properties[key]

	try:
		# NOT GENERIC ENOUGH!!!
		data = HTT_generator(legs, jets, ints, floats)
		# print(next(data))
		
		# writing to file
		binary_template.wpisz(data)

	except:
	    print "[ERROR] GENERATION FAILED! MESSAGE:", sys.exc_info()[0]
	    raise

def read_binary():
	dataset = binary_template.wczytaj_dataset()
	return dataset
	# BATCH_SIZE=2
	# zbachowany=dataset.shuffle(1000).repeat().batch(BATCH_SIZE)
	# iterator = zbachowany.make_one_shot_iterator()
	# f,l=iterator.get_next()
	# print("tu dzial")
	# with tf.Session() as sess:
	#     for i in range(10):
	#         print(sess.run([f,l]))

if __name__ == "__main__":
    save2binary()
    read_binary()


	
	
