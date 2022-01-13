#!/usr/bin/env python3
import regex as re
import pandas as pd
import numpy as np



def CE_MismatchCapture(PrimerSeq, OptimialBindingDict):

	# Complementary Dictionary - can expand to contain
	# special characters as well.
	CompDict = {'A':['T'],
				'T':['A'],
				'C':['G'],
				'G':['C']}

	# Only generate classification string if mismatches present in anchor
	# region. (i.e. 0 = mismatch and 1 = successful binding)

	print('Binding string')
	print(OptimialBindingDict['MaxRelBindngSeq'])

	if 0 in OptimialBindingDict['MaxRelBindngSeq'][-4:]:

		# Convert primer-template anchor region into classification string
		# ????-???? for mismatch characterisation.
		PrimerAnchor = PrimerSeq[-4:]
		TemplateAnchor = OptimialBindingDict['FinalSubSeq3to5'][-4:]

		PStr = ''
		TStr = ''

		for i in range(4):

			# Check if complementary
			if PrimerAnchor[i] in CompDict[TemplateAnchor[i]]:
				PStr+='?'
				TStr+='?'

			else:
				PStr+=PrimerAnchor[i]
				TStr+=TemplateAnchor[i]

		# Combine to form classification string
		MMClass = 'MM_Anchor:{}-{}'.format(PStr,TStr)

		return MMClass

	else:
		return False
