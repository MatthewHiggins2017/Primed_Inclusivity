#!/usr/bin/env python3


# Function to capture length of primer from the 3' terminus binding
# desired target site, terminating with 3 or more adjacent mismatches.
def CE_MinLenCapture(PrimerSeq, OptimialBindingDict):

	# Extract binding array, and strip possible leading mismatches from
	# 5' terminus. Then convert back into list and invert so now 3'-5'
	BindingArray = list(''.join([str(x) for x in OptimialBindingDict['MaxRelBindngSeq']]).lstrip('0'))[::-1]


	# Determine min-length from 3' terminus.
	MismatchCount = 0
	LastValidLength = 0
	for i in range(len(BindingArray)):
		if BindingArray[i] == '1':
			LastValidLength = i
			MismatchCount = 0

		else:
			MismatchCount +=1
			if MismatchCount == 3:
				break

	# Have to +1 as python zero indexed.
	MinLenEntry = 'MinLen:{}'.format(str(LastValidLength+1))


	return MinLenEntry
