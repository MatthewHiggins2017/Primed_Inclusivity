#!/usr/bin/env python3
import regex as re
import pandas as pd
import numpy as np


def CE_TotalBindingCount(PrimerSeq,
					     OptimialBindingDict):

	# Count Matches
	Matches = OptimialBindingDict['MaxRelBindngSeq'].count(1)

	return 'MinLen:{}'.format(str(Matches))



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
