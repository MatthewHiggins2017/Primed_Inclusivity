#!/usr/bin/env python3
import numpy as np
import regex as re
import sympy as sp
import pandas as pd


# Function to determine the probability according to mismatch anchor.
def MM_Anchor(ClassificationKeys, MMAnchorGuide, SC):

	# Extract anchor classification from list
	# which should be of length 1.
	AnchorString = ClassificationKeys[0].split(':')[-1]

	# For MMProbContents, the keys for the MM_Anchor dictionary described in the
	# json should match the different Anchor strings. Therefore, first need to
	# check if the Anchor String is present in Guide
	if AnchorString in list(MMAnchorGuide.keys()):

		# Return the selected dictionary which contains values for P,E,OT etc.
		return MMAnchorGuide[AnchorString]


	# Convert to positional classification and extract the corresponding
	# metrics from the hardcoded list below.
	else:
		PositionalOutputDict = {'T':{'P':0.876,
									 'E':19.366,
									 'OT':17.715},

								'T1N':{'P':0.811,
									 'E':40.144,
									 'OT':34.172},

								'T2N':{'P':0.824,
									 'E':32.810,
									 'OT':28.891},

								'T3N':{'P':0.781,
									 'E':32.100,
									 'OT':27.352},

								 '1N':{'P':0.965,
 									 'E':9.666,
 									 'OT':1.758},

 								 '2N':{'P':0.960,
 									  'E':13.057,
 									  'OT':-0.707},

								 '3N':{'P':0.940,
 									  'E':15.326,
 									  'OT':1.072}
								}

		# Boolean indicates if mismsatch is present or not.
		PositionalDict = {'3N':False,
						  '2N':False,
						  '1N':False,
						  'T':False}

		# For each position determine if mismatch is present or not.
		for P, i in zip(PositionalDict.keys(),range(len(PositionalDict))):
			if AnchorString[i]!='?':
				PositionalDict[P]=True

		# Generate mismatch positional string (arranged in correct way)
		PositionalString = ''.join([a for a,b in PositionalDict.items() if b == True][::-1])


		return PositionalOutputDict[PositionalString]
