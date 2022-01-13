import json
import numpy as np
from pkgs.Output import *
from collections import defaultdict

################################################################################

# Given classification string determine the probability
# basically build a logic system which will encode.
def OutputEngine(ClassificationString,
				 ToolParameters,
			     SetConditions,
				 ):

	# Probability list
	FinalProbList = []

	# Split Classification string into list of items
	ClassificationList = ClassificationString.split(',')
	ClasIndList = [c.split(':')[0] for c in ClassificationList]

	# Load in Output Engine Guide
	with open(ToolParameters.OutputGuide) as Prob_json_file:
		ProbJson = json.load(Prob_json_file)


	# Define dictionary which contains the relevant output variables
	OutputStore = defaultdict(list)

	# Extract Output Functions from JSON
	# Loop through key indicators
	for Indicator, Contents in ProbJson["Output Engine"]['Output Functions'].items():

		# Check if Indicator actually in the classification string
		if Indicator in ClasIndList:

			# Subset classification list to contain only those with
			# variables with indicator of interest.
			ClassificationOfInterest = [C for C in ClassificationList if Indicator in C]

			# Run function named according to the indicator and pass through the
			# relevant values. This will return a dictionary with the derived values
			# corresponding to each output value e.g {'E':100,'P':0.8}. Make sure
			# all variables return as dictionary. Fix this as hard rule.

			OutputDict = eval(Indicator + '({},{},{})'.format(ClassificationOfInterest,
															 Contents,
															 SetConditions))

			# Add Outputs to the final OutputStore Store
			for OI, OV in OutputDict.items():
				OutputStore[OI] += [OV]


	# Generate the final output values.
	# and store as dictionary. If an output is missing then
	# assign as reference value according to output json

	FinalOutput = {}
	for OI in list(ToolParameters.OutputVar.keys()):
		if OI in OutputStore.keys():

			FinalOutput[OI] = CombineOutput(OutputStore[OI], ToolParameters.OutputVar[OI]['Binding_Combination'])
		else:
			FinalOutput[OI] = ToolParameters.OutputVar[OI]['Reference_Value']


	return FinalOutput

################################################################################


def CombineOutput(SL,KeyPhrase):


	CSL = [float(l) for l in SL]


	Output = None
	KeyPhrase = KeyPhrase.strip().upper()

	if KeyPhrase == "MIN":
		Output = np.min(CSL)

	elif KeyPhrase == "MAX":
		Output = np.max(CSL)

	elif KeyPhrase== "MEAN":
		Output = np.mean(CSL)

	elif KeyPhrase== "PRODUCT":
		Output = np.product(CSL)

	elif KeyPhrase== "SUM":
		Output = np.sum(CSL)

	return Output
