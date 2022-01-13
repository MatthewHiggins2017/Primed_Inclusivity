#!/usr/bin/env python3
import numpy as np
import regex as re
import sympy as sp
import pandas as pd




def Thermo(ThermoString ,BFormula, SetConditions):
	# Extract Thermo Values from Classification String.
	TD = {}
	TD['Gibbs'], TD['Enthalpy'], TD['Entropy'],  PhosphatesPresent, ReliabilityScore = [float(px) for px in ThermoString[0].split(':')[1:]]


	# Performed Set-Specific Salt Adjustments Adjustments
	TDU, MT, BindingFrac = ThermoSaltCorrection(TD,
												PhosphatesPresent,
												SetConditions,
												RunningTm=SetConditions['Annealing_Temperature'])



	F = BFormula['Formula']
	X = sp.symbols('X')
	BFExpression = sp.sympify(F)

	Prob = BFExpression.subs('X', BindingFrac)

	print('BindingFrac',BindingFrac,'Prob',Prob)
	return {'P':Prob}






def ThermoSaltCorrection(ThermoDict,
						 PhosphatesPresent,
						 ReagentConcs,
						 RunningTm=False):

	# Calculate Salt Molarity for penalisation and turn user input which is
	# in mM to M.
	if ReagentConcs['Divalent_Cation_Concentration_(mM)'] > ReagentConcs['Total_dNTP_Concentration_(mM)']:
		SaltConc = (ReagentConcs['Monovalent_Cation_Concentration_(mM)'] + (120*((ReagentConcs['Divalent_Cation_Concentration_(mM)']-ReagentConcs['Total_dNTP_Concentration_(mM)'])**0.5))) * 10**-3
	else:
		SaltConc = ReagentConcs['Monovalent_Cation_Concentration_(mM)'] * 10**-3


	# Add in Salt Adjustment
	ThermoDict['Adjusted Gibbs'] =  ThermoDict['Gibbs'] + (-0.114 * (PhosphatesPresent/2) * np.log(SaltConc))
	ThermoDict['Adjusted Entropy'] = ThermoDict['Entropy'] + (0.368 * (PhosphatesPresent/2) * np.log(SaltConc))


	# Define Primer Concentration total converting
	# from uM to Moles.
	CT = ReagentConcs['Primer_Concentration_(uM)']*10**-6


	# Determine Melting Temperature
	R=1.9872
	MeltingTemp = (ThermoDict['Enthalpy']*1000)/(ThermoDict['Adjusted Entropy']+(R*np.log((CT)/4)))-273.15


	# Determine Percentage Binding if RunningTm is not equal to false.
	# The section below is taken from the Juypter Notebook.
	if RunningTm != False:

		H = ThermoDict['Enthalpy']
		S = ThermoDict['Adjusted Entropy']
		T = RunningTm

		KD =  1 / np.exp(((H * 1000)/ (R * (T + 273.15))) - (S/R))


		### Formula for determinine Ratio of split based on temperature
		Top = (KD * CT) - np.sqrt(((2 * KD * CT) + 1)) + 1
		Bottom = (KD * CT)
		BindingFraction = Top/Bottom

	else:
		BindingFraction = np.nan

	print('MeltingTemp',MeltingTemp, 'BindingFraction', BindingFraction)
	return (ThermoDict, MeltingTemp, BindingFraction)
