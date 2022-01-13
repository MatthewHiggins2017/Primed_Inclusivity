import os
import sys
import ast
import json
import sqlite3
import itertools
import subprocess
import regex as re
import numpy as np
import sympy as sp
import pandas as pd
from ipfn import ipfn
import multiprocessing
from collections import ChainMap, defaultdict
from .ClassificationEngine import *
from .OutputEngine import *


# Create SQL Database
def create_sql_connection(db_file):
	conn = None
	try:
		conn = sqlite3.connect(db_file)
		return conn
	except Error as e:
		print(e)
	return conn


# Directory Check
def CheckAndMake(Path):
	if os.path.exists(Path):
		return True
	else:
		os.mkdir(Path)
		return True


# File Check
def FilePresent(Path,Kill):
	if os.path.exists(Path) == False:
		if Kill == True:
			print("Error: Incorrect path for: {}".format(Path))
			sys.exit()
		else:
			return False



# VCF Sample List Extraction
def VCFSampleExtract(VCFPath):
	ExtractSamples = 'bcftools query -l {}'.format(VCFPath)
	SamplesListQuery = subprocess.run(ExtractSamples,shell=True,capture_output=True)
	SampleList = SamplesListQuery.stdout.decode("utf-8").strip().split('\n')

	return SampleList


# Obtain Complement Sequence
def getComplement(seq,
				 reverse=False,
				 regex=False,
				 rule='N2N'):

	seqComp = ""
	for base in seq.upper():
		base = base.upper()
		if base == "A":
			seqComp += "T"
		elif base == "T":
			seqComp += "A"
		elif base == "C":
			seqComp += "G"
		elif base == "G":
			seqComp += "C"
		elif base == 'U':
			seqComp += "A"
		elif base == 'R':
			seqComp += "Y"
		elif base == 'Y':
			seqComp += "R"
		elif base == 'S':
			seqComp += "S"
		elif base == 'W':
			seqComp += "W"
		elif base == 'K':
			seqComp += "M"
		elif base == 'M':
			seqComp += "K"
		elif base == '.':
			if regex == True:
				seqComp += '.'
			else:
				seqComp += 'N'
		elif base == "N":
			if rule == 'N2-':
				seqComp += '-'
			elif rule == 'N2N':
				seqComp += "N"
		elif base == "-":
			seqComp += "N"

		# If Character not any of above assign it as N
		else:
			seqComp += "N"

	if reverse:
		return(seqComp[::-1])
	else:
		return(seqComp)



# Determine the Alt Binding Seequences Extraction.
def BindingSiteExtraction(PrimerID, ToolParameters):

	# Establish New Connection and Cursor Objects, ready for row extraction
	TempConn = create_sql_connection(ToolParameters.DatabasePath)
	TempConn.row_factory = sqlite3.Row
	TempCur = TempConn.cursor()

	# Extract Row of Interest Containing Primer Information
	TempCur.execute("SELECT * FROM Primer_Table WHERE  PrimerID=:ID",
				   {"ID": PrimerID})
	RowOfI = TempCur.fetchone()


	# Now extend primer binding search equivalent to 1 primer length
	# downstream and 1 primer length upstream. This will really help
	# in finding the optimal binding string.

	BSStartIndex =	RowOfI['BindingSite'] - len(RowOfI['PrimerSeq'])
	BSEndIndex = RowOfI['BindingSite'] + (len(RowOfI['PrimerSeq'])*2)

	if BSStartIndex < 0:
		BSStartIndex = 0


	# Subset and extract VCF key columns to generate alt binding sequence. - DOUBLE CHECK +1 is necessary.
	GenGTVCFCommand = "bcftools query -f '%CHROM:%POS:%REF:%ALT \n' -r {0}:{1}-{2} {3} -S {4}".format(RowOfI['BindingChromo'],
																								BSStartIndex + 1,
																								BSEndIndex + 1,
																								ToolParameters.VCF,
																								ToolParameters.SamplesFile)

	GenGTVCF = subprocess.run(GenGTVCFCommand,shell=True,capture_output=True)
	GenGTVCFString = GenGTVCF.stdout.decode("utf-8")


	GenGTVCFDF = pd.DataFrame([R.split(':') for R in GenGTVCFString.split('\n')],
							   columns=['CHROM','POS','REF','ALT']).dropna().drop_duplicates()

	# Set Index as POS Column and subtract one to make inline with
	# python 0 indexing
	GenGTVCFDF.index = [int(PI)-1 for PI in GenGTVCFDF['POS'].tolist()]


	# Check to if the first variants is less than the binding site start
	# site as if so set this as the start index.
	if GenGTVCFDF.index.tolist()[0] < BSStartIndex:
		BSStartIndex = GenGTVCFDF.index.tolist()[0]



	# Extract Reference + sense 5'-3' sequence including 1x primer length
	# up and downstream.
	RefWideBindingSeq = NucleotideFinder(ToolParameters.Target,
										 RowOfI['BindingChromo'],
										 BSStartIndex,
										 1+BSEndIndex-BSStartIndex)


	# Determine Sample Clusters Around GT Strings
	SampleClusters = ExtractSampleClusters(ToolParameters,
										   RowOfI['BindingChromo'],
										   BSStartIndex + 1,
										   BSEndIndex + 1)



	# Define Sample Clusters where the GT string is completely the WT i.e.
	# all individual variable sites correspond to the reference i.e 0:0:0.
	AllRefGTString = ':'.join(['0'] * len(list(SampleClusters.keys())[0].split(':')))


	# If AllRefGTString not in SampleClusters.keys() add it in
	# so perfect bindinng will be captured and classified accordingly.
	if AllRefGTString not in SampleClusters.keys():
		SampleClusters[AllRefGTString] = []


	# Extract and sort all GT strings. By sorting ensure the WT with perfect binding
	# goes first.
	SortGTSL = sorted(list(SampleClusters.keys()))

	KVR = list(RowOfI.keys())

	CLROWI = {KVR[x]:y for x,y in enumerate(RowOfI)}

	GenABInput = zip(SortGTSL,
					 [ToolParameters]*len(SampleClusters),
					 [GenGTVCFDF]*len(SampleClusters),
					 [CLROWI]*len(SampleClusters),
					 [BSStartIndex]*len(SampleClusters),
					 [RefWideBindingSeq]*len(SampleClusters)
					 )


	with multiprocessing.Pool(ToolParameters.Threads) as Pool:
		AltBSMetrics = Pool.starmap(CharacteriseAltBindingSeq,GenABInput)


	# Identify GT Strings which result in the same binding information and then
	# I can index according to unique alternative binding sites.

	GTrev = defaultdict(list)
	for pda in AltBSMetrics:
		GTS = list(pda.keys())[0]
		InfoDict = str(list(pda.values())[0])
		GTrev[InfoDict]+=[GTS]


	# Need to ensure Perfect Binding is first, so that 0 index
	# always corresponds to perfect binding.
	GTRevInfo = list(GTrev.keys())
	HOP = []
	for GI in range(len(GTRevInfo)):
		HOP = GTRevInfo[GI]
		BA = ast.literal_eval(HOP)['Binding Array']
		# Confirm if all Binding Array is all 1 i.e. perfect binding
		# then make sure its at start of list
		if set(BA) == {1}:
			GTRevInfo.pop(GI)
			GTRevInfo = [HOP] + GTRevInfo
			break

	# Now we have the characteristics for each GT string and what samples contain
	# each GT string. I need to assign index from 1 onways to alt binding seq as
	# remember 0 = perfect binding.
	GTIndexList = range(0,len(GTrev))


	# Create Sample default dict to store information
	# regarding sample classification.
	SampleGTID = defaultdict(list)


	# Update SQL Alt Binding Table according to each alternative binding seq.
	for GTIndex in GTIndexList:

		GTCharacterisation = ast.literal_eval(GTRevInfo[GTIndex])


		AltBindingSeqUpdate = '''INSERT INTO Binding_Seq_Table(PrimerID, BindingSeqIndex,
								 BindingSeq, BindingArray, ClassificationString) VALUES(?,?,?,?,?)'''

		AltBindingSeqvalues = (PrimerID,
							   GTIndex,
							   GTCharacterisation['Alt Binding Seq'],
							   ','.join([str(GTB) for GTB in GTCharacterisation['Binding Array']]),
							   GTCharacterisation['Binding Classification']
							   )

		TempConn.execute(AltBindingSeqUpdate,AltBindingSeqvalues)
		TempConn.commit()


		# Loop through samples in the relevant cluster for the GT String.
		# and update the default dictionary

		for GTStringOfI in GTrev[GTRevInfo[GTIndex]]:

			for SOI in SampleClusters[GTStringOfI]:
				SampleGTID[SOI]+=[str(GTIndex)]


	# Update the Sample vs Primer SQL Table to identify the variant present
	# for each sample. Later down the line if the cell is empty i.e. null
	# it means that it equals the reference and should be included as 0.
	for Sample, BindingList in SampleGTID.items():

		BindingString = ','.join(BindingList)

		# SQL Updated
		SampleVsPrimerUpdate = 	"UPDATE Sample_vs_Primer_Table SET {} = ? WHERE Samples = ?".format(PrimerID)

		TempConn.execute(SampleVsPrimerUpdate,(BindingString,Sample.strip()))
		TempConn.commit()

	# Return True once all primer information is derived
	# and associated sql tables updated.
	return True



# Given a VCF determine the cluster of samples
# present according to GT present. VALIDATE and TEST!!!!!
def ExtractSampleClusters(ToolParameters,
						  ChromoOfIntr,
						  StartIndex,
						  EndIndex):


	# Subset VCF, remove header, and convert into pandas dataframe.
	SubsetVCFCommand = "bcftools view -H -r {0}:{1}-{2} {3} -S {4}".format(ChromoOfIntr,
																	StartIndex,
																	EndIndex,
																	ToolParameters.VCF,
																	ToolParameters.SamplesFile) # Parse Class to function


	SubsetVCFCommandQuery = subprocess.run(SubsetVCFCommand,shell=True,capture_output=True)
	SubsetVCFString = SubsetVCFCommandQuery.stdout.decode("utf-8")
	SubsetVCFDF = pd.DataFrame([LP.split('\t') for LP in SubsetVCFString.split('\n')]).dropna()


	# Extract VCF header and add to subset dataframe.
	ExtractVCFHeader = 'bcftools view -h {} -S {}'.format(ToolParameters.VCF,ToolParameters.SamplesFile)

	ExtractVCFHeaderQuery = subprocess.run(ExtractVCFHeader,shell=True,capture_output=True)
	ExtractVCFHeaderString = ExtractVCFHeaderQuery.stdout.decode("utf-8")
	HLL = [LL for LL in ExtractVCFHeaderString.split('\n') if '#CHROM' in LL][0]
	VCFHeaderIndex = HLL.split('\t')


	SubsetVCFDF.columns = VCFHeaderIndex

	# Subset VCF to contain Sample Genotypes
	SampleGTs = SubsetVCFDF.iloc[:,9:].T

	# This dictionary will store the unique GT strings and Sample IDs
	UniqueGTStrings = defaultdict(list)

	# Loop through samples and create unique GTs
	for SID in SampleGTs.index.tolist():

		# Determine the GT Strings present for each sample
		# accounting for phased chromosomes
		SampleSpecGT = defaultdict(list)

		# Loop through variant positions.
		for varpos in SampleGTs.loc[SID,:]:


			# EDITING !!!! - Errors here which need fixing
			SplitGT = varpos.replace('.','0').replace('/','|').split(':')[0].split('|')


			for ChromoIndex in range(len(SplitGT)):
				SampleSpecGT[ChromoIndex+1] += [SplitGT[ChromoIndex]]

		# Update the UniqueGTStrings default dict
		for Chromo, GTSL in SampleSpecGT.items():
			GTString = ':'.join(GTSL)
			UniqueGTStrings[GTString]+= [SID]

	troll = UniqueGTStrings.keys()

	# Now we have obtained all unique GT strings across the population
	# and can assign chromosome clusters. However the only thing to take care of
	# is to ensure no duplicate sample names for GT Strings
	CleanUniquGTDict = {gtstring:sorted(list(set(samplelist))) for gtstring, samplelist in UniqueGTStrings.items()}

	return CleanUniquGTDict




# VALIDATE AND TEST!!!!!
def CharacteriseAltBindingSeq(GTString,
							  ToolParameters,
							  GuideVCF,
							  PrimerRowInfo,
							  BSStartIndex,
							  RefWideBindingSeq):


	# Split GT String Back Into list of integers
	GTList = [int(GT) for GT in GTString.split(':')]
	GTTracker = 0


	# Establish where Alt Binding Seq will be written.
	AltBindingSeq = ''


	# Loop through positions and update the AltBindingSeq accordingly
	NS = 0

	# Excess to add on if deletions.
	EAD = 0

	# Determine Starting Position according to start minus 1x primer length
	#BSStartIndex = PrimerRowInfo['BindingSite'] - len(PrimerRowInfo['PrimerSeq'])


	# While steps of NS not equal to RefWideBindingSeq as this is to account
	# for excess sequence up and downstream to allow full overlap assessment
	while NS <= len(RefWideBindingSeq)-1:


		# Determine relavative binding position.
		BSOI = BSStartIndex + NS

		# Check if variant present at site. Remember python is 0 indexed
		# whereas VCF is 1 index therefore position 70 in VCF corresponds
		# to 69 in the python sequence

		if BSOI in GuideVCF.index:

			# Extract Key Information
			VarInfo = GuideVCF.loc[BSOI,:].to_dict()

			# Extract Relevant Variant and add to Binding String
			VarOptions = [VarInfo['REF']] + [av for av in VarInfo['ALT'].split(',')]

			VarOI = VarOptions[GTList[GTTracker]].strip()
			AltBindingSeq += VarOI

			# Update Step based on length of reference.
			NS += len(VarInfo['REF'])

			# Update Excess based on differences
			if len(VarInfo['REF']) >= len(VarOI):
				EAD += len(VarInfo['REF']) - len(VarOI)

			# Update GTTracker
			GTTracker += 1


		# If no variant present at site update the AltBindingSeq using same
		# nucleotides as the primer of interest.
		else:
				AltBindingSeq += RefWideBindingSeq[NS]
				NS+=1

			# Else going to have to look to fasta to fill in the remaining positions
			# this is key when there is a deletion. This will always finish off the
			# alt binding sequence.

	# If Excess Present due to deletions fill it.
	if EAD != 0:


		# Determine nucleotides needed to fill the excess.
		LastNucls = NucleotideFinder(ToolParameters.Target,
									 PrimerRowInfo['BindingChromo'],
									 BSStartIndex + NS, # This should equal end of seq
									 EAD)


		# Update AltBindingSeq to complete it
		AltBindingSeq += LastNucls


	# Now AltBindingSeq is created move onto determine optimal binding and
	# classify. However before doing so make sure to correct sequences according
	# to sense etc to make sure I am optimising the binding etc.

	# Alt Biding seq has been created to be 5' to 3' + Sense
	# therefore if primer is + sense need to complement the
	# alt binding seq to make 3'-5' direction. However if primer is - sense will reverse alt
	# binding seq to read in 3'-5; direction WRITE OUT TO MAKE SURE THIS MAKES SENSE !!!!

	if PrimerRowInfo['Sense'] == '+':
		FinalAltBindingSeq = getComplement(AltBindingSeq,reverse=False)
	else:
		FinalAltBindingSeq = AltBindingSeq[::-1]


	# FOR THIS FUNCTION ALWAY HAVE THE PRIMER FIRST IN 5'-3' on top AND THEN THE
	# ALT BINDING SEQ SECOND 3'-5'. RETURN A DICT BASED OBJECT AS WILL BE EASIER
	# TO MANIPULATE. POSSIBLE KEEP SIMPLE IN TERMS OF MAXIMUM BINDING SITES. ALSO
	# POTENTIALLY ADD IN SHIFT EXPANSION TO ALT BINDING SEQ.
	OptimalBindingDict = FindOptimalBinding(PrimerRowInfo['PrimerSeq'],
											FinalAltBindingSeq,
											ToolParameters)


	# Generate Binding Site Specific Classification String.
	ClassificationString = BindingSiteClassify(PrimerRowInfo['PrimerSeq'],
											   OptimalBindingDict)



	# Format output ready to return and to be filled into
	# the alt binding sequences dataframe.

	return { GTString: { 'Alt Binding Seq':OptimalBindingDict['FinalSubSeq3to5'],
						 'Binding Array':OptimalBindingDict['MaxRelBindngSeq'],
						 'Binding Classification':ClassificationString}}



def FindOptimalBinding(Seq5to3,
					   Seq3to5,
					   OBParam):


	# Start with both zero indexed and shift across recording maximum
	# matches. Assumes Seq3to5 is longer than Seq5to3
	Seq5to3 = Seq5to3.upper()
	Seq3to5 = Seq3to5.upper()


	OptimalBindingDict = {'Weighting':-1000}

	# Expand to include IUPAC Codes
	ComplementaryBinding = {'A':['T'],
							'T':['A'],
							'C':['G'],
							'G':['C']}


	# Subset Seq3to5 according to each step and length of Seq5to3 then
	# run comparison.
	SI = 0
	for SubSeq3to5 in [Seq3to5[i:i+len(Seq5to3)] for i in range(len(Seq3to5)-len(Seq5to3))]:

		# Make faster later -- clean up.
		RelBindingSeq = []
		MatchCounter = 0
		for Pos in range(len(Seq5to3)):
			if Seq5to3[Pos] in ComplementaryBinding[SubSeq3to5[Pos]]:
				MatchCounter += 1
				RelBindingSeq.append(1)
			else:
				RelBindingSeq.append(0)


		# Using Rel Binding Seq apply rules and weightings to
		# determine the optimal binding.
		Weighting = RelBindingSeq.count(1) * OBParam.ComplementaryScalar

		# Determine length of rel binding sequence.
		RBSLen = len(RelBindingSeq)

		# Loop through positions.
		for rp in range(len(RelBindingSeq)):
			AdjacentMismatch = False

			# Determine if position of interest is mismatched.
			if (RelBindingSeq[rp] == 0):

				# Determine if adjacent position is mismatched.
				if (rp + 1 <= RBSLen - 1) and (RelBindingSeq[rp+1] == 0):
					AdjacentMismatch = True

				# If adjacent mismatch true and
				# mismatch within anchor region.
				if (AdjacentMismatch == True) and (rp >= RBSLen-OBParam.AnchorRegion):
					Weighting += -1  * OBParam.AnchorRegionMismatchScaler * \
									   OBParam.AnchorAdjacentMismatchScalar

				# If adjacent mismatch true and
				# outside anchor region.
				elif (AdjacentMismatch == True) and (rp < RBSLen-OBParam.AnchorRegion):
					Weighting += -1  * OBParam.MismatchScaler * \
									   OBParam.AdjacentMismatchScalar

				# If adjacent mismatch false and
				# single mismatch outside anchor region
				elif (AdjacentMismatch == False) and (rp < RBSLen-OBParam.AnchorRegion):
					Weighting += -1  * OBParam.MismatchScaler

				# If adjacent mismatch false and
				# single mismatch inside anchor region
				elif (AdjacentMismatch == False) and (rp > RBSLen-OBParam.AnchorRegion):
					Weighting += -1  * OBParam.AnchorRegionMismatchScaler

		if Weighting > OptimalBindingDict['Weighting']:
			OptimalBindingDict['Weighting'] = float("%.5f" % Weighting)
			OptimalBindingDict['MaxMatches'] = MatchCounter
			OptimalBindingDict['ShiftIndex'] = SI
			OptimalBindingDict['FinalSubSeq3to5'] = SubSeq3to5
			OptimalBindingDict['MaxRelBindngSeq'] = RelBindingSeq

		# Update step
		SI +=1


	return OptimalBindingDict



# VALIDATE AND TEST!!!!!!
def NucleotideFinder(FastaPath,
					 ChromoID,
					 StartPos,
					 NucleoidesNeeded):

	ChromoFound = False
	PosCounter = 0

	ExtractFound = False
	Extract = ''

	with open(FastaPath, 'rb') as f:
		for line in f:
			line = line.decode("utf-8").replace('\n','')
			if '>{}'.format(ChromoID) in line:
				ChromoFound = True
				continue


			# Incase extract goes over multiple lines need to place this assessment
			# before extract exception so it will move onto next line if necessary.
			if (len(Extract) != NucleoidesNeeded) and (ExtractFound == True):
				Extract+=line[:NucleoidesNeeded-len(Extract)]


			if (ChromoFound == True) and (ExtractFound == False):
				if len(line) + PosCounter < StartPos:
					PosCounter += len(line)
				else:
					LineRelPos = StartPos - PosCounter
					Extract += line[LineRelPos:LineRelPos+NucleoidesNeeded]
					ExtractFound = True


			# Return Extract
			if (len(Extract) == NucleoidesNeeded) and (ExtractFound == True):
				return Extract.upper()




# TEST and VALIDATE
# Extract Set specific probability and extract key information.
def DeriveSetMetrics(SetID, ToolParameters):

	TempConn = create_sql_connection(ToolParameters.DatabasePath)
	TempConn.row_factory = sqlite3.Row
	TempCur = TempConn.cursor()


	# Extract Row of Interest Containing SetID.
	TempCur.execute("SELECT * FROM Sets_Table WHERE  SetID=:ID",
				   {"ID": SetID})
	SetRowOfI = TempCur.fetchone()


	# This should return dictionary object as I have created a
	# row factory.


	# Extract Primers of Interest
	RelevantPrimers = SetRowOfI['Primers_Included'].split(':')


	# Extract Set Conditions, and remove hard-coded fixed columns
	SetConditionsDict = {k:SetRowOfI[k] for k in SetRowOfI.keys() if k not in ['index','SetID','PrimerIncluded']}


	# Create new set-specific SQL table by subsetting the Sample vs Primer
	# table.
	SetSpecDataBaseId = 'Set_{}_S_v_P'.format(SetRowOfI[1])

	#CreateSetSpecTable = '''CREATE TABLE {} AS SELECT * FROM Sample_vs_Primer_Table WHERE PrimerID IN ({})'''.format(SetSpecDataBaseId,
	#																												','.join(['?']*len(RelevantPrimers)))


	# SQL UPDATED
	CreateSetSpecTable = '''CREATE TABLE {} AS SELECT {} FROM Sample_vs_Primer_Table'''.format(SetSpecDataBaseId,','.join(RelevantPrimers+['Samples']))


	TempCur.execute(CreateSetSpecTable)

	# For the Set specific table update all the values to reflect the probabilies.
	# For the primer of interest derive set specific probabilities.
	# and update the Set Specific Sample vs Primer test.
	for PID in RelevantPrimers:

		DeriveSetSpecProb(PID,
						  SetSpecDataBaseId,
						  ToolParameters,
						  SetConditionsDict)



	# Now derive the summary sample values for the given primer set by
	# averaging all relevant primer information.

	for OVI in list(ToolParameters.OutputVar.keys()):

		NCS = 'Sample_{}'.format(OVI)

		# Add new columns to table, ready to hold summary values.
		AddNewColumns = '''ALTER TABLE {} ADD {} VARCHAR(100);'''.format(SetSpecDataBaseId,NCS)

		TempConn.execute(AddNewColumns)
		TempConn.commit()


	for RSample in ToolParameters.SamplesList:

		# Extract Sample Probs across Set Primers
		COLextraction = 'SELECT {} FROM {} WHERE Samples = (?)'.format(','.join(RelevantPrimers),SetSpecDataBaseId)
		TempCur.execute(COLextraction, [RSample])
		SampleOutStrings = TempCur.fetchall()

		# Extract the values according to string into list form.
		SampleOutDict = defaultdict(list)
		for SOS in SampleOutStrings[0]:
			for OC in SOS.split(','):
				OI, OV = OC.split(':')
				SampleOutDict[OI]+=[OV]


		# Determine mean output value for each sample
		# for each output variable.
		for OVI in list(ToolParameters.OutputVar.keys()):

			SSOV = CombineOutput(SampleOutDict[OVI], ToolParameters.OutputVar[OVI]['Set_Combination'])

			SampleVsPrimerUpdate = 	"UPDATE {} SET Sample_{} = ? WHERE Samples = ?".format(SetSpecDataBaseId, OVI)
			TempConn.execute(SampleVsPrimerUpdate,[SSOV, RSample])
			TempConn.commit()

	return True




def DeriveSetSpecProb(PrimerID,
					  SetTableOfInterest,
					  ToolParameters,
					  SetConditions):


	# Extract row of interest from the table
	TempConn = create_sql_connection(ToolParameters.DatabasePath)
	TempConn.row_factory = sqlite3.Row
	TempCur = TempConn.cursor()

	# Extract PrimerID Row
	TempCur.execute("SELECT {} FROM {}".format(PrimerID, SetTableOfInterest))
	RowOfI = TempCur.fetchall()
	RowOfI = [x[0] for x in RowOfI]

	#Extract Sample IDs
	TempCur.execute("SELECT Samples FROM {}".format(SetTableOfInterest))
	ColumnNames = TempCur.fetchall()
	ColumnNames = [x[0] for x in ColumnNames]

	# For each Sample determine the alt binding site indexes present, derive the
	# the corresponding probabilities and then keep the top probabity. To save time
	# create a dictionary
	PrimerSetSpecProbDict = defaultdict(float)

	for SI in range(len(ColumnNames)):
		Sample = ColumnNames[SI]
		ProCd = RowOfI[SI]

		if Sample not in ['index','PrimerID']:

			#SampleSpecValues = []

			PFS = defaultdict(list)

			# If binding site index for a given primer sample pair has not been
			# assigned this means it is the reference and so the binding site is
			# 0.
			if ProCd is not None:
				if len(ProCd) != 0:
					AltBindingSiteIndexes = ProCd.split(',')
				else:
					AltBindingSiteIndexes = ['0']

			else:
				AltBindingSiteIndexes = ['0']

			for ABSI in AltBindingSiteIndexes:

				if ABSI not in PrimerSetSpecProbDict.keys():

					# Look up binding index in Binding_Seq_Table
					# the generate corresponding probability. For now
					# limit this to the Classification string.

					CSExtr = 'SELECT ClassificationString FROM Binding_Seq_Table WHERE BindingSeqIndex = (?) AND PrimerID = (?)'
					TempCur.execute(CSExtr,[ABSI,PrimerID])
					CSString = TempCur.fetchone()[0]


					OutDict = OutputEngine(CSString,
										   ToolParameters,
										   SetConditions)



					PrimerSetSpecProbDict[ABSI] = OutDict
				else:
					OutDict = PrimerSetSpecProbDict[ABSI]

				# Add to values to necessary variable entries in dictionary.
				for OI,OV in OutDict.items():
					PFS[OI] += [OV]



			FinalStringList = []
			for FOI,FOV in PFS.items():

				CFV = CombineOutput(FOV,ToolParameters.OutputVar[FOI]['Sample_Combination'])

				FinalStringList.append('{}:{}'.format(FOI,str(CFV)))

			FinalString = ','.join(FinalStringList)

			# Update the Set Specific Sample vs Primer table. To contain the final string
			SampleVsPrimerUpdate = 	"UPDATE {} SET {} = ? WHERE Samples = ?".format(SetTableOfInterest,
																					 PrimerID)


			TempConn.execute(SampleVsPrimerUpdate,(FinalString,Sample))
			TempConn.commit()

	return True




def GenerateRepresentationTable(RepresentationPath):

	'''
	Pre-requisite:
		- The order in which the classifiers are defined with respective
		  marginal distributions has too be the same order for the classifiers
		  of the samples.
	'''

	with open(RepresentationPath) as Rep_json_file:
		Repjson = json.load(Rep_json_file)


	RepTags = Repjson['Representation']['RepTags']
	SampleClassification = Repjson['Representation']['SampleTags']

	AllCVars = [list(cv.keys()) for cv in RepTags.values()]
	AllCVcombos = [':'.join(CVC) for CVC in itertools.product(*AllCVars)]
	CVCIndexes = list(range(len(AllCVcombos)))
	CVCRefDict = dict(zip(AllCVcombos, CVCIndexes))


	SamplesPerRefCombo = defaultdict(list)


	for SampleID, CVTagsDict in SampleClassification.items():
		SampleSpecCVStr = ':'.join(CVTagsDict.values())
		RelRepID = CVCRefDict[SampleSpecCVStr]
		SamplesPerRefCombo[RelRepID] += [SampleID]


	# Now have all representation combinations and also corresponding samples
	# I can now convert this into a DataFrame

	RDCols =['Rep_Index'] + list(RepTags.keys()) + ['Samples', 'total']
	RepDF = pd.DataFrame(columns=RDCols)

	for ComboString, RepIndex in CVCRefDict.items():

		RelSamp = ':'.join(SamplesPerRefCombo[RepIndex])
		RelSampCount = len(SamplesPerRefCombo[RepIndex])
		NR = [RepIndex] + ComboString.split(':') + [RelSamp, RelSampCount]
		RepDF.loc[len(RepDF)] = NR


	# Now I have RepDF with all useful information and it is in a similar format
	# to the https://pypi.org/project/ipfn/ pandas example. However i need intelligent
	# way to do dimension finding even when I dont know the relevant agregates and
	# dimensions.

	AggDict = {}
	Dimensions = []

	for CV in list(RepTags.keys()):
		Dimensions.append([CV])
		TempLA = RepDF.groupby(CV)['total'].sum()
		# Extract the actual marginal distributions and update TempLA
		for TLAK in TempLA.keys():
			TempLA[TLAK] = RepTags[CV][TLAK]

		AggDict[CV] = TempLA


	Aggregates = [AggDict[CVT] for CVT in list(RepTags.keys())]

	# Archive the unadjusted total count values. Then run the
	# Iterative proportional fitting algorithm to adjust totals.

	RepDF['total_unadjusted'] = RepDF['total']
	IPF = ipfn.ipfn(RepDF, Aggregates, Dimensions)
	UpdRepDF = IPF.iteration()



	# Next steps will be to add the sets information as columns and subsequently
	# two additional rows. Where the Rep_Index actually equals a string such
	# as unadjusted overall probability, adjusted overall probability for one table
	# then unadjusted perfect count and adjusted perfect count for the other
	# table. Is it worth including this in this function though or as an
	# additional functions.

	return UpdRepDF













def IdentifyBindingSitesV2(OligoRegexString,
						   reference,
						   Overlap):


	ChromoDict = {}
	BindingSites = []

	# Keep Track of Index and Chromo
	active_chromo = 0
	Index = 0
	CapturedOverlap = ''


	with open(reference) as r:
		for line_idx, line in enumerate(r):
			line = line.strip().replace(' ','').replace('\n','').upper()
			if not line: # Check line is not empty
				continue

			# Update Chromo Dict to keep track of potential
			# binding sites.
			if line.startswith(">"):
				active_chromo +=1
				ChromoDict[line[1:]] = active_chromo

				# Reset Index so Chromosome specific
				Index = 0
				CapturedOverlap = ''

				continue

			# Search For Possible Binding Sites
			else:

				# Merge line with prior carryover.
				WorkingLine = CapturedOverlap + line

				# Iterate over Both Sense of Working line as now keeping
				# Oligo Regex string fixed so need to reverse complement the Working Line,
				# when necessary.

				for sense in ['F','R']:
					if sense == 'R':
						TempWorkingLine = getComplement(WorkingLine,True)
					else:
						TempWorkingLine = WorkingLine

					for OligoIndex in  [m.start(0) for m in re.finditer(OligoRegexString, TempWorkingLine, overlapped=True)]:

						# Extract the respecitve Oligo binding index accounting
						# for sense differences.
						if sense == 'F':
							AdjTrueIndex = OligoIndex + Index
						else:
							AdjTrueIndex = (len(TempWorkingLine)-OligoIndex) + Index - len(OligoRegexString)

						# Append to the Binding Site list.
						BindingSites.append((active_chromo,sense,AdjTrueIndex))

				# Update captured overlap and adjust the index accordingly
				#CapturedOverlap = line[-Overlap:]
				CapturedOverlap = WorkingLine[-Overlap:]
				Index += len(WorkingLine)-len(CapturedOverlap)


	return (ChromoDict, list(set(BindingSites)))



def DeriveCombinedOutput(OutputVar,
						 ToolParameters):

	# Extract row of interest from the table
	TempConn = create_sql_connection(ToolParameters.DatabasePath)
	TempConn.row_factory = sqlite3.Row
	TempCur = TempConn.cursor()
	TempCur.execute("SELECT * FROM {0}_Representation_Table WHERE Rep_Index NOT IN (?, ?)".format(OutputVar),
				   ['Unadjusted_{}'.format(OutputVar),'Adjusted_{}'.format(OutputVar)])
	RL = TempCur.fetchall()

	# Loop through each Rep_Index present
	for ROI in RL:
		SampleIDs = ROI['Samples'].split(':')

		# Derive and loop through set lists
		for SetID in ToolParameters.SetIDList:

			# Conncet to Set specific tables and then extract all Sample_Combined_Probability
			# or other dependent variable values where Samples equal to values in the
			# list defined above. Then take mean
			EP = 'SELECT Sample_{} FROM Set_{}_S_v_P WHERE Samples IN ({})'.format(OutputVar,
																		   	  	   SetID,
																		   	  	   ','.join(['?']*len(SampleIDs)))

			TempCur.execute(EP,SampleIDs)
			EPAll = TempCur.fetchall()


			# Extract all Sample specific variable values, ready to calculate the
			# mean value.
			RSSP = [float(PR['Sample_{}'.format(OutputVar)]) for PR in EPAll]



			RSSPUC = "UPDATE {}_Representation_Table SET `{}` = ? WHERE Rep_Index == ?".format(OutputVar,SetID)


			TempCur.execute(RSSPUC,[np.mean(RSSP),ROI['Rep_Index']])
			TempConn.commit()




	# Now the relevant represenation table is updated with all corresponding values.
	# I can then perform set based calculations. Again this needs to be done on a
	# set specific level. Also derive unadjusted and adjusted values

	# Extract Representation totals and convert to array
	ERT = 'SELECT total FROM {}_Representation_Table WHERE Rep_Index NOT IN (?, ?)'.format(OutputVar)
	TempCur.execute(ERT, ['Unadjusted_{}'.format(OutputVar),'Adjusted_{}'.format(OutputVar)])


	RepArray = np.array([RR['total'] for RR in TempCur.fetchall()])

	# Dirty Fix = Polish Later - This is so unadjusted values match that if the
	# representation file was not included.
	if ToolParameters.RepresentationGuide != '':
		DirtyERT = 'SELECT total_unadjusted FROM {}_Representation_Table WHERE Rep_Index NOT IN (?, ?)'.format(OutputVar)
		TempCur.execute(DirtyERT, ['Unadjusted_{}'.format(OutputVar),'Adjusted_{}'.format(OutputVar)])
		UnadjRepArray = np.array([RR['total_unadjusted'] for RR in TempCur.fetchall()])
	else:
		UnadjRepArray = RepArray


	for SetID in ToolParameters.SetIDList:

		SSID = 'SELECT `{}` FROM {}_Representation_Table WHERE Rep_Index NOT IN (?, ?)'.format(SetID,OutputVar)
		TempCur.execute(SSID, ['Unadjusted_{}'.format(OutputVar),'Adjusted_{}'.format(OutputVar)])

		TPArray = np.array([float(RR[0]) for RR in TempCur.fetchall()])

		for UC in ['Unadjusted','Adjusted']:

			if UC == 'Unadjusted':
				V = np.sum(UnadjRepArray * TPArray)/np.sum(UnadjRepArray)
			else:
				V = np.sum(RepArray * TPArray)/np.sum(RepArray)


			URDFC = "UPDATE {}_Representation_Table SET `{}` = ? WHERE Rep_Index = ?".format(OutputVar,SetID)

			TempCur.execute(URDFC,(V,'{}_{}'.format(UC,OutputVar)))
			TempConn.commit()


			# Update the final sets table as well, so user can look at this for
			# simply summary output without all the information required for the
			# derivation.
			USTC = 'UPDATE Sets_Table SET `{}` = ? WHERE SetID = ?'.format('{}_{}'.format(UC,OutputVar))

			TempCur.execute(USTC,(V,SetID))
			TempConn.commit()


	return True


def ExtractOutputVar(Outputjson):
	with open(Outputjson) as OJ:
		OJE = json.load(OJ)
	return OJE['Output Engine']['Output Variables']
