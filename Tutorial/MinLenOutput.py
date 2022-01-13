


def MinLen(ClassificationKeys, MinLenGuide, SC):

	# Extract MinLen classification from list
	# which should be of length 1.
	MinLen = int(ClassificationKeys[0].split(':')[-1])

	# Acording to the MinLength guide determine if the value
	# obtained falls below or above the threshold for a successful
	# reaction and set output variable accordingly
	ReactionPass = {'RS':1}
	if MinLen < MinLenGuide['Threshold']:
		ReactionPass['RS'] = 0

	return ReactionPass
