# Purpose = To test core functions to ensure perform as expected.

from General.Core import *
from shutil import copyfile
import filecmp



# Switch off tests
'''

#==============================================================================#
################################################################################
# NucleotideFinder Tests
################################################################################

FastaPath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test.fa'
ChromoID = 'Test_Driven_Dev'
StartPos = 5
NucleoidesNeeded = 10

NFTest1 = NucleotideFinder(FastaPath,
                           ChromoID,
                           StartPos,
                           NucleoidesNeeded)

if NFTest1 == 'CTCTGGTTAG':
    print('Pass')
else:
    print('Fail NucleotideFinder() Test 1')
    sys.exit()

################################################################################

FastaPath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test2.fa'
ChromoID = 'Chromo2'
StartPos = 5
NucleoidesNeeded = 10

NFTest2 = NucleotideFinder(FastaPath,
                           ChromoID,
                           StartPos,
                           NucleoidesNeeded)

if NFTest2 == 'TAGACAAATA':
    print('Pass')
else:
    print('Fail NucleotideFinder() Test 2')
    sys.exit()



################################################################################


FastaPath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test2.fa'
ChromoID = 'Chromo2'
StartPos = 5
NucleoidesNeeded = 50

NFTest2 = NucleotideFinder(FastaPath,
                           ChromoID,
                           StartPos,
                           NucleoidesNeeded)

if NFTest2 == 'TAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAG':
    print('Pass')
else:
    print('Fail NucleotideFinder() Test 3')
    sys.exit()



################################################################################
#==============================================================================#





#==============================================================================#
################################################################################
# FindOptimalBinding() Tests
################################################################################


Seq5to3 = 'AATTGGCC'
Seq3to5 = 'TTTTTTTTTTTTTTTAACCGGAAAAAAAAAAAAA'

DesiredDict = {'MaxMatches': 8, 'ShiftIndex': 13, 'FinalSubSeq3to5': 'TTAACCGG', 'MaxRelBindngSeq': [1, 1, 1, 1, 1, 1, 1, 1]}

if DesiredDict == FindOptimalBinding(Seq5to3,
				                     Seq3to5):
    print('Pass')
else:
    print('Fail FindOptimalBinding() Test 1')
    sys.exit()

################################################################################


Seq5to3 = 'AATCGGCC'
Seq3to5 = 'TTTTTTTTTTTTTTTAACCGGAAAAAAAAAAAAA'

DesiredDict = {'MaxMatches': 7, 'ShiftIndex': 13, 'FinalSubSeq3to5': 'TTAACCGG', 'MaxRelBindngSeq': [1, 1, 1, 0, 1, 1, 1, 1]}

if DesiredDict == FindOptimalBinding(Seq5to3,
				                     Seq3to5):
    print('Pass')
else:
    print('Fail FindOptimalBinding() Test 2')
    sys.exit()

################################################################################
#==============================================================================#




#==============================================================================#
################################################################################
# BindingSiteClassify() Tests
################################################################################


PrimerSeq = 'AATCGGCC'
OptimalBindingDict =  {'MaxMatches': 7, 'ShiftIndex': 13, 'FinalSubSeq3to5': 'TTAACCGG', 'MaxRelBindngSeq': [1, 1, 1, 0, 1, 1, 1, 1]}

if 'MM_Count:1,4:C-A' == BindingSiteClassify(PrimerSeq,
					                          OptimalBindingDict):
    print('Pass')
else:
    print('Fail BindingSiteClassify() Test 1')
    sys.exit()


################################################################################

PrimerSeq = 'AATCGGCA'
OptimalBindingDict =  {'MaxMatches': 6, 'ShiftIndex': 13, 'FinalSubSeq3to5': 'TTAACCGG', 'MaxRelBindngSeq': [1, 1, 1, 0, 1, 1, 1, 0]}

if 'MM_Count:2,0:A-G,4:C-A' == BindingSiteClassify(PrimerSeq,
					                                OptimalBindingDict):
    print('Pass')
else:
    print('Fail BindingSiteClassify() Test 2')
    sys.exit()

################################################################################
#==============================================================================#





#==============================================================================#
################################################################################
# ExtractSampleClusters() Tests
################################################################################

VCFPath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Haploid_Test.vcf.gz'
ChromoOfIntr = 'Test_Driven_Dev'
StartIndex = 60
EndIndex = 90


ExpectedOutputDict = {'1:1': ['Sample_1', 'Sample_3'], '0:0': ['Sample_2', 'Sample_4']}


if ExpectedOutputDict == ExtractSampleClusters(VCFPath,
					                           ChromoOfIntr,
                    					        StartIndex,
                    					        EndIndex):
    print('Pass')
else:
    print('Fail ExtractSampleClusters() Test 1')
    sys.exit()

################################################################################
# Add another advanced tets with more advanced VCF
################################################################################
#==============================================================================#





#==============================================================================#
################################################################################
# AltBindingSiteExtraction() Tests
################################################################################

class ToolParam:
    DatabasePath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test.db'
    Target = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test.fa'
    VCF = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Haploid_Test.vcf.gz'
    Threads = 1


copyfile('/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Original_Test.db',
        '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test.db')

PrimerID = 'FP_1'

if __name__ == '__main__':
    AltBindingSiteExtraction(PrimerID,ToolParam)


# Check output test.db gainst Test_val_1.db to make sure identical

# Sample Clustering
# Sample 1 and 3 = 1:1
# Sample 2 and 4 = 0:0 and therefore should be dropped as matches reference

# Region which is extracted.
# GCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAG

# Sample 1 and 2 should be:
#GCTAACTAGGGAACCCACTGCTTAAGCTTCTAAAGCTTGCCTTGAGTGCTTCAAGTA
#GCTAACTAGGGAACCCACTGCTTAAGCTTCTAAAGCTTGCCTTGAGTGCTTCAAGTA


# Sample 2 and 4 should match the reference.
#GCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAG
#GCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAG

# To validate check if database matches and at the start
# copy database and update path so the testing does not
# break the test database and can continue to test repeatedly.


################################################################################
# Add another advanced test
################################################################################
#==============================================================================#





#==============================================================================#
################################################################################
# ProbabilityEngine()  + MM() Test
################################################################################

# Should find positions 0 and 1 and then use regex probabilities of 0.8 to extract
# therefore the combined output prob would be 0.8 * 0.8 = 0.64

ClassificationString = 'MM_Count:2,MM_Pos:0:G-G,MM_Pos:1:G-G'
ProbJsonGuide = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Probability_Guide_val_1.json'

vv = ProbabilityEngine(ClassificationString,
				        ProbJsonGuide)

x ="%.2f" % vv

if x == '0.64':
    print('Pass')
else:
    print('Fail ProbabilityEngine() Test 1')
    sys.exit()


################################################################################


ClassificationString = 'MM_Count:2,MM_Pos:0:A-C,MM_Pos:1:G-G'
ProbJsonGuide = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Probability_Guide_val_1.json'

vv = ProbabilityEngine(ClassificationString,
				        ProbJsonGuide)



x ="%.3f" % vv

if x == '0.008':
    print('Pass')
else:
    print('Fail ProbabilityEngine() Test 2')
    sys.exit()






################################################################################
#==============================================================================#


'''

#==============================================================================#
################################################################################
# DeriveSetMetrics()
################################################################################

input('DeriveSetMetrics Test')


# With the starting database Follow_Up_Test the RP perfectly binds to all samples
# 1-4, whereas the FP only binding to samples 2 and 4 perfectly. For samples 1 and 3
# the FP has the following classification string MM_Count:7,0:C-A,1:G-A,2:A-G,3:A-C,5:T-T,7:A-A,10:C-A.
# When using the Probability_Guide_val_1.json this means that both full under the regex
# mismsatch classifier and so successful amplificaiton from the primer should be
# classified as 0.64, for the FP for samples 1 and 3. Also as the probability of
# the RP is 1, then the probability for samples 1-4 respectively should be 0.64, 1, 0.64,
# 0.64. Whereby the average probability is then 0.82 and the percentage of perfect
# binding is 50%.




copyfile('/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Follow_Up_Test.db',
        '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Set_Metric_Test.db')


class ToolParam:
    DatabasePath = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Set_Metric_Test.db'
    Target = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Test.fa'
    Threads = 1
    ProbJsonGuide = '/Users/matthewhiggins/Documents/PhD/Programming/Github/PrimedInclusivity/Validation/TDD/Probability_Guide_val_1.json'


SetID ='1'

DeriveSetMetrics(SetID,
                 ToolParam)




################################################################################
#==============================================================================#
