#!/usr/bin/env python3
import numpy as np
import regex as re
import sympy as sp
import pandas as pd


# Function to determine the probability from classification keys
# provided according to mismatch.
def MM_Pos(ClassificationKeys, MMProbContents, SC):

	OutDict = {}

	ProbabilityList = []
	EfficiencyList = []
	OnsetList = []

	# Convert classification keys into dictionary
	# i.e. 'MM:0:A-C' becomes {'0':['A','C']}
	CleanClassKeys = {}
	for K in ClassificationKeys:
		Pos, Nucs = K.replace('MM_Pos:','').split(':')
		CleanClassKeys[Pos] = Nucs.split('-')


	for TermPos, NucInvovled in CleanClassKeys.items():
		# Check position exists in keys
		if TermPos in MMProbContents.keys():
			# Loop through options in position till match is found
			for PosOptions in MMProbContents[TermPos]:
				if ((bool(re.match('N{}'.format(PosOptions['Primer_Nucleotide']), 'N{}'.format(NucInvovled[0])))==True)&
				   (bool(re.match('N{}'.format(PosOptions['Template_Nucleotide']),'N{}'.format(NucInvovled[1])))==True)):
					ProbabilityList.append(PosOptions['P'])
					EfficiencyList.append(PosOptions['E'])
					OnsetList.append(PosOptions['OT'])
					break

	# Once all probability keys have been found take the product of the list
	if len(ProbabilityList) != 0:
		OutDict['P'] = np.prod(ProbabilityList)
	if len(EfficiencyList) != 0:
		OutDict['E'] = np.sum(EfficiencyList)
	if len(EfficiencyList) != 0:
		OutDict['OT'] = np.sum(OnsetList)


	return OutDict


def MinLen(MLString,MLGuideDict,SC):
	print('\n\n MLString',MLString)
	OutDict = {}
	if int(MLString[0].split(':')[-1].strip()) <= 20:
		 OutDict['P'] = MLGuideDict['Fail']
	return OutDict
