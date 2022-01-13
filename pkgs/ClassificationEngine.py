import pkgs.Classification

################################################################################

def BindingSiteClassify(PrimerSeq,
						OptimialBindingDict):

	FunctionsInPlay = [F for F in dir(pkgs.Classification) if 'CE_' == F[:3]]

	ClassificationList = []

	for CF in FunctionsInPlay:
		CSString = eval('pkgs.Classification.'+CF + '("{}",{})'.format(PrimerSeq,OptimialBindingDict))
		if CSString != False:
			ClassificationList.append(CSString)

	ClassificationString = ','.join(ClassificationList)

	return ClassificationString


################################################################################
