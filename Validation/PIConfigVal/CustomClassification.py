#!/usr/bin/env python3
import regex as re
import pandas as pd
import numpy as np


def CE_TotalMismatchCount(PrimerSeq,
					   OptimialBindingDict):

	# Count mismatches
	Mismatches = OptimialBindingDict['MaxRelBindngSeq'].count(0)

	return 'MM_Count:{}'.format(str(Mismatches))



def CE_MismatchCapture(PrimerSeq, OptimialBindingDict):

	# Capture mismatch positions walking from 3' end
	# where by 0 = terminal nucleotide in 3' end and
	# 1 = 1n etc etc. In the string include the nucleotides
	# involved in the mismatch.
	TempCl = []

	for ttp in range(len(OptimialBindingDict['MaxRelBindngSeq'])):
		if OptimialBindingDict['MaxRelBindngSeq'][::-1][ttp] == 0:
			TempCl.append('MM_Pos:{}:{}-{}'.format(str(ttp),
												   PrimerSeq[::-1][ttp],
												   OptimialBindingDict['FinalSubSeq3to5'][::-1][ttp]))


	return ','.join(TempCl)



# Modify to derive actual thermodynmamic values.
def CE_ThermoClassify(PrimerSeq,
				   OptimialBindingDict):


	# Determine if john santalucia binding is satisfied and obtain the
	# forward and reverse oligos.

	TP, ModOpBD = ThermoCompatibleRegion(PrimerSeq,
										 OptimialBindingDict)

	if TP == True:
		ThermoDict = Thermodynamics(ModOpBD['ThermoBindingSeqs'][0],
									ModOpBD['ThermoBindingSeqs'][1])


		return 'Thermo:{}:{}:{}:{}:{}'.format(ThermoDict['Gibbs'],
											  ThermoDict['Enthalpy'],
											  ThermoDict['Entropy'],
											  ModOpBD['ThermoPhosphatePresent'],
											  ModOpBD['ThermoReliability']
												)

	# If no compatible Binding then simply return as false so no
	# thermo tag is includes in the classification string.
	else:
		return False


def ThermoCompatibleRegion(PrimerSeq,
						  OptimalBindingDict):

	# Inverse RelBindingSeq so going 3' to 5' direction
	InvRelBinding = OptimalBindingDict['MaxRelBindngSeq'][::-1]

	# Check if T and 1N are not both mismatched and if they are then
	# no Thermo Binding can be defined as anchor region is not set.
	if (InvRelBinding[0]==0) and (InvRelBinding[1]==0):
		return (False, OptimalBindingDict)


	# Head forward and find binding according to rules. Can have terminal
	# mismatches each end and no adjacent mismatches in middle. Expanding
	# down from full length.
	for CL in range(len(InvRelBinding)):

		PosBinding = InvRelBinding[0:len(InvRelBinding)-CL]

		# Check no two 0 are adjacent using regex and that the region is at
		# least big enough according to user for thermodynamic calculation.
		RegexPassCheck = True

		# https://stackoverflow.com/questions/60962752/pattern-to-detect-repeating-characters-not-working-regex
		# The regex search below will find strings with
		# adjacent 0 (mismatch) characters.
		if (re.search(r'(0)\1', ''.join([str(M) for M in PosBinding]))):
			RegexPassCheck = False


		MinSize = 5

		if (RegexPassCheck == True) and (CL <= len(InvRelBinding)-MinSize):

			OptimalBindingDict['ThermoBindingDuplex']=True
			OptimalBindingDict['ThermoBindingSeqs'] = (PrimerSeq[CL:],
													   OptimalBindingDict['FinalSubSeq3to5'][CL:])

			OptimalBindingDict['ThermoPhosphatePresent']  = (len(PrimerSeq[CL:])*2)-2


			# Determine % of binding outside of thermodynamic limit according to
			# the size of the fragment not include. I.e. need to ensure reliability
			# is not compromised. i.e. if double adjacent mismatch then full binding this
			# will throw out the binding potential of primers.

			ExcludedRegion = InvRelBinding[len(InvRelBinding)-CL:]

			ExcludedRegionLen = len(ExcludedRegion)

			if ExcludedRegionLen != 0:

				ExcludedRegionMatches = ExcludedRegion.count(1)


				OptimalBindingDict['ThermoReliability'] = (ExcludedRegionMatches /  ExcludedRegionLen) * \
														  (ExcludedRegionLen / len(InvRelBinding))

			else:
				OptimalBindingDict['ThermoReliability'] = 1


			return (True, OptimalBindingDict)


	# If the condition is never satisfied i.e. too short
	# or always adjacent mismatches.
	return (False, OptimalBindingDict)




def Thermodynamics(PlusSense,MinusSense):

	# Define ThermoDF
	ThermoDF = pd.DataFrame({'Pair One': {0: 'AA',
	  1: 'AT',
	  2: 'TA',
	  3: 'CA',
	  4: 'GT',
	  5: 'CT',
	  6: 'GA',
	  7: 'CG',
	  8: 'GC',
	  9: 'GG',
	  10: 'AA',
	  11: 'CA',
	  12: 'GA',
	  13: 'TA',
	  14: 'AC',
	  15: 'CC',
	  16: 'GC',
	  17: 'TC',
	  18: 'AG',
	  19: 'CG',
	  20: 'GG',
	  21: 'TG',
	  22: 'AT',
	  23: 'CT',
	  24: 'GT',
	  25: 'TT',
	  26: 'AG',
	  27: 'AT',
	  28: 'CG',
	  29: 'CT',
	  30: 'GG',
	  31: 'GG',
	  32: 'GT',
	  33: 'GT',
	  34: 'TG',
	  35: 'TG',
	  36: 'TT',
	  37: 'AA',
	  38: 'AC',
	  39: 'CA',
	  40: 'CC',
	  41: 'GA',
	  42: 'GC',
	  43: 'TA',
	  44: 'TC',
	  45: 'AA',
	  46: 'AG',
	  47: 'CA',
	  48: 'CG',
	  49: 'GA',
	  50: 'GG',
	  51: 'TA',
	  52: 'TG',
	  53: 'AC',
	  54: 'AT',
	  55: 'CC',
	  56: 'CT',
	  57: 'GC',
	  58: 'GT',
	  59: 'TC',
	  60: 'TT'},
							 'Pair Two': {0: 'TT',
							  1: 'TA',
							  2: 'AT',
							  3: 'GT',
							  4: 'CA',
							  5: 'GA',
							  6: 'CT',
							  7: 'GC',
							  8: 'CG',
							  9: 'CC',
							  10: 'TA',
							  11: 'GA',
							  12: 'CA',
							  13: 'AA',
							  14: 'TC',
							  15: 'GC',
							  16: 'CC',
							  17: 'AC',
							  18: 'TG',
							  19: 'GG',
							  20: 'CG',
							  21: 'AG',
							  22: 'TT',
							  23: 'GT',
							  24: 'CT',
							  25: 'AT',
							  26: 'TT',
							  27: 'TG',
							  28: 'GT',
							  29: 'GG',
							  30: 'CT',
							  31: 'TT',
							  32: 'CG',
							  33: 'TG',
							  34: 'AT',
							  35: 'GT',
							  36: 'AG',
							  37: 'TC',
							  38: 'TA',
							  39: 'GC',
							  40: 'GA',
							  41: 'CC',
							  42: 'CA',
							  43: 'AC',
							  44: 'AA',
							  45: 'TG',
							  46: 'TA',
							  47: 'GG',
							  48: 'GA',
							  49: 'CG',
							  50: 'CA',
							  51: 'AG',
							  52: 'AA',
							  53: 'TT',
							  54: 'TC',
							  55: 'GT',
							  56: 'GC',
							  57: 'CT',
							  58: 'CC',
							  59: 'AT',
							  60: 'AC'},
							 'Delta H': {0: -7.6,
							  1: -7.2,
							  2: -7.2,
							  3: -8.5,
							  4: -8.4,
							  5: -7.8,
							  6: -8.2,
							  7: -10.6,
							  8: -9.8,
							  9: -8.0,
							  10: 1.2,
							  11: -0.9,
							  12: -2.9,
							  13: 4.7,
							  14: 0.0,
							  15: -1.5,
							  16: 3.6,
							  17: 6.1,
							  18: -3.1,
							  19: -4.9,
							  20: -6.0,
							  21: 1.6,
							  22: -2.7,
							  23: -5.0,
							  24: -2.2,
							  25: 0.2,
							  26: 1.0,
							  27: -2.5,
							  28: -4.1,
							  29: -2.8,
							  30: 3.3,
							  31: 5.8,
							  32: -4.4,
							  33: 4.1,
							  34: -0.1,
							  35: -1.4,
							  36: -1.3,
							  37: 2.3,
							  38: 5.3,
							  39: 1.9,
							  40: 0.6,
							  41: 5.2,
							  42: -0.7,
							  43: 3.4,
							  44: 7.6,
							  45: -0.6,
							  46: -0.7,
							  47: -0.7,
							  48: -4.0,
							  49: -0.6,
							  50: 0.5,
							  51: 0.7,
							  52: 3.0,
							  53: 0.7,
							  54: -1.2,
							  55: -0.8,
							  56: -1.5,
							  57: 2.3,
							  58: 5.2,
							  59: 1.2,
							  60: 1.0},
							 'Delta S': {0: -21.3,
							  1: -20.4,
							  2: -21.3,
							  3: -22.7,
							  4: -22.4,
							  5: -21.0,
							  6: -22.2,
							  7: -27.2,
							  8: -24.4,
							  9: -19.9,
							  10: 1.7,
							  11: -4.2,
							  12: -9.8,
							  13: 12.9,
							  14: -4.4,
							  15: -7.2,
							  16: 8.9,
							  17: 16.4,
							  18: -9.5,
							  19: -15.3,
							  20: -15.8,
							  21: 3.6,
							  22: -10.8,
							  23: -15.8,
							  24: -8.4,
							  25: -1.5,
							  26: 0.9,
							  27: -8.3,
							  28: -11.7,
							  29: -8.0,
							  30: 10.4,
							  31: 16.3,
							  32: -12.3,
							  33: 9.5,
							  34: -1.7,
							  35: -6.2,
							  36: -5.3,
							  37: 4.6,
							  38: 14.6,
							  39: 3.7,
							  40: -0.6,
							  41: 14.2,
							  42: -3.8,
							  43: 8.0,
							  44: 20.2,
							  45: -2.3,
							  46: -2.3,
							  47: -2.3,
							  48: -13.2,
							  49: -1.0,
							  50: 3.2,
							  51: 0.7,
							  52: 7.4,
							  53: 0.2,
							  54: -6.2,
							  55: -4.5,
							  56: -6.1,
							  57: 5.4,
							  58: 13.5,
							  59: 0.7,
							  60: 0.7},
							 'Gibbs free energy': {0: -1.0,
							  1: -0.88,
							  2: -0.58,
							  3: -1.45,
							  4: -1.44,
							  5: -1.28,
							  6: -1.3,
							  7: -2.17,
							  8: -2.24,
							  9: -1.84,
							  10: 0.61,
							  11: 0.43,
							  12: 0.17,
							  13: 0.69,
							  14: 1.33,
							  15: 0.7,
							  16: 0.79,
							  17: 1.05,
							  18: -0.13,
							  19: -0.11,
							  20: -1.11,
							  21: 0.44,
							  22: 0.69,
							  23: -0.12,
							  24: 0.45,
							  25: 0.68,
							  26: 0.71,
							  27: 0.07,
							  28: -0.47,
							  29: -0.32,
							  30: 0.08,
							  31: 0.74,
							  32: -0.59,
							  33: 1.15,
							  34: 0.43,
							  35: 0.52,
							  36: 0.34,
							  37: 0.88,
							  38: 0.77,
							  39: 0.75,
							  40: 0.79,
							  41: 0.81,
							  42: 0.47,
							  43: 0.92,
							  44: 1.33,
							  45: 0.14,
							  46: 0.02,
							  47: 0.03,
							  48: 0.11,
							  49: -0.25,
							  50: -0.52,
							  51: 0.42,
							  52: 0.74,
							  53: 0.64,
							  54: 0.73,
							  55: 0.62,
							  56: 0.4,
							  57: 0.62,
							  58: 0.98,
							  59: 0.97,
							  60: 0.75},
							 'Resource': {0: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  1: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  2: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  3: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  4: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  5: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  6: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  7: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  8: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  9: 'Table 1 - doi: 10.1146/annurev.biophys.32.110601.141800',
							  10: '10.1021/bi9825091',
							  11: '10.1021/bi9825091',
							  12: '10.1021/bi9825091',
							  13: '10.1021/bi9825091',
							  14: '10.1021/bi9825091',
							  15: '10.1021/bi9825091',
							  16: '10.1021/bi9825091',
							  17: '10.1021/bi9825091',
							  18: '10.1021/bi9825091',
							  19: '10.1021/bi9825091',
							  20: '10.1021/bi9825091',
							  21: '10.1021/bi9825091',
							  22: '10.1021/bi9825091',
							  23: '10.1021/bi9825091',
							  24: '10.1021/bi9825091',
							  25: '10.1021/bi9825091',
							  26: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  27: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  28: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  29: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  30: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  31: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  32: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  33: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  34: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  35: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  36: 'Table 5 - https://doi.org/10.1021/bi962590c',
							  37: 'Table 4 -  DOI: 10.1021/bi9803729',
							  38: 'Table 4 -  DOI: 10.1021/bi9803729',
							  39: 'Table 4 -  DOI: 10.1021/bi9803729',
							  40: 'Table 4 -  DOI: 10.1021/bi9803729',
							  41: 'Table 4 -  DOI: 10.1021/bi9803729',
							  42: 'Table 4 -  DOI: 10.1021/bi9803729',
							  43: 'Table 4 -  DOI: 10.1021/bi9803729',
							  44: 'Table 4 -  DOI: 10.1021/bi9803729',
							  45: 'Table 4 - DOI: 10.1021/bi9724873',
							  46: 'Table 4 - DOI: 10.1021/bi9724873',
							  47: 'Table 4 - DOI: 10.1021/bi9724873',
							  48: 'Table 4 - DOI: 10.1021/bi9724873',
							  49: 'Table 4 - DOI: 10.1021/bi9724873',
							  50: 'Table 4 - DOI: 10.1021/bi9724873',
							  51: 'Table 4 - DOI: 10.1021/bi9724873',
							  52: 'Table 4 - DOI: 10.1021/bi9724873',
							  53: 'Table 4 - 10.1093/nar/26.11.2694',
							  54: 'Table 4 - 10.1093/nar/26.11.2694',
							  55: 'Table 4 - 10.1093/nar/26.11.2694',
							  56: 'Table 4 - 10.1093/nar/26.11.2694',
							  57: 'Table 4 - 10.1093/nar/26.11.2694',
							  58: 'Table 4 - 10.1093/nar/26.11.2694',
							  59: 'Table 4 - 10.1093/nar/26.11.2694',
							  60: 'Table 4 - 10.1093/nar/26.11.2694'}
							  })

	# Define Dictionary
	OutputDict = {'Gibbs':1.96,
				  'Entropy':-5.7,
				  'Enthalpy':0.2,
				  'Gibbs5P':1.96,
				  'Gibbs3P':1.96}


	GibbsFreeEnergyInput = list(zip([PlusSense[i:i+2] for i in list(range(len(PlusSense)-1))],
									[MinusSense[i:i+2] for i in list(range(len(MinusSense)-1))]))


	# Iterate through Nearest-Neighbour paired nucleotides.
	MissingPair = False
	for NNpair in GibbsFreeEnergyInput:
		try:
			RowOfInterestIndex = ThermoDF[(ThermoDF['Pair One']==NNpair[0])&(ThermoDF['Pair Two']==NNpair[1])].index.tolist()[0]
		except:
			try:
				RowOfInterestIndex = ThermoDF[(ThermoDF['Pair One']==NNpair[1][::-1])&(ThermoDF['Pair Two']==NNpair[0][::-1])].index.tolist()[0]
			except:
				print(NNpair)
				MissingPair = True

		# If Nearest-Neighbour pair valid update thermodynamic scores.
		if MissingPair == False:
			OutputDict['Gibbs'] += ThermoDF.loc[RowOfInterestIndex,'Gibbs free energy']
			OutputDict['Enthalpy'] += ThermoDF.loc[RowOfInterestIndex,'Delta H']
			OutputDict['Entropy'] += ThermoDF.loc[RowOfInterestIndex,'Delta S']


	# If NN pair not found, return function early.
	# and generate necessary error message.
	if MissingPair == True:
		return ('Fail NN Pair not Found')


	# Derive terminal gibbs free energy values
	FivePrime = True
	for TerminalPairs in [GibbsFreeEnergyInput[:6],
						  GibbsFreeEnergyInput[-6:]]:

		for PLNN in TerminalPairs:
			try:
				II = ThermoDF[(ThermoDF['Pair One']==PLNN[0])&
							  (ThermoDF['Pair Two']==PLNN[1])].index.tolist()[0]
			except:
				II = ThermoDF[(ThermoDF['Pair One']==PLNN[1][::-1])&
							  (ThermoDF['Pair Two']==PLNN[0][::-1])].index.tolist()[0]


			if FivePrime == True:
				OutputDict['Gibbs5P'] += ThermoDF.loc[II,'Gibbs free energy']
			else:
				OutputDict['Gibbs3P'] += ThermoDF.loc[II,'Gibbs free energy']

		FivePrime = False


	# Apply AT penality if necessary.
	FivePrime = True
	for terminal in [''.join([PlusSense[0],MinusSense[0]]),
					 ''.join([PlusSense[-1],MinusSense[-1]])]:

		if terminal in ['AT','TA']:


			OutputDict['Gibbs'] += 0.05
			OutputDict['Enthalpy'] += 2.2
			OutputDict['Entropy'] +=  6.9

			if FivePrime == True:
				OutputDict['Gibbs5P'] += 0.05
				FivePrime = False
			else:
				OutputDict['Gibbs3P'] += 0.05

	return OutputDict
