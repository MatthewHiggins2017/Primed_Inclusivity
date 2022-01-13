#!/usr/bin/env python3

def CE_PerfectBinding(PrimerSeq,
				   OptimialBindingDict):

	# Determine if perfect binding. I.e. if a mismatch
    # is present then perfect binding is false.
    if 0 in OptimialBindingDict['MaxRelBindngSeq']:
        return 'PB:0'
    else:
        return 'PB:1'
