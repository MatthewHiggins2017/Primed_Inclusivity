import os
import json
import subprocess
import pandas as pd
import filecmp
from pandas._testing import assert_frame_equal
# Extract location of the PrimedInclusivity Software
PackagePath = os.path.dirname(os.path.abspath(__file__))


#################################################################################
# PIConfig Validation
#################################################################################

print("\n\nPIConfig Validation\n\n")
input('Enter')


# Create Custom Module Json File
CustomJson = {"PIConfig":{
                         "Classification":['{}/pkgs/Classification/DefaultClassification.py'.format(PackagePath),
                                           '{}/Validation/PIConfigVal/CustomClassification.py'.format(PackagePath)],

                         "Output":['{}/pkgs/Output/DefaultOutput.py'.format(PackagePath),
                                   '{}/Validation/PIConfigVal/CustomOutput.py'.format(PackagePath)],
                            }}

CustomJsonPath = '{}/Validation/PIConfigVal/PIConfigVal.json'.format(PackagePath)

with open(CustomJsonPath, 'w') as CJP:
    json.dump(CustomJson, CJP)
CJP.close()



# Run PIConfig Command Custom
PIConfigCustomCommand = 'python PIConfig -CF {}'.format(CustomJsonPath)
subprocess.run(PIConfigCustomCommand,shell=True)



# Check Init Files match and pass if they do.

if (filecmp.cmp('{}/pkgs/Output/__init__.py'.format(PackagePath),
                  '{}/Validation/PIConfigVal/OC__init__.py'.format(PackagePath))) == True:

    print('Output Modules Pass')

if (filecmp.cmp('{}/pkgs/Classification/__init__.py'.format(PackagePath),
                  '{}/Validation/PIConfigVal/CC__init__.py'.format(PackagePath))) == True:

    print('Classification Modules Pass')




#Run PIConfig Command to return to default state.
PIDefaultCustomCommand = 'python PIConfig -D'
subprocess.run(PIDefaultCustomCommand,shell=True)


#################################################################################
# ValQuick Automated Run
#################################################################################
print("\n\nValidation Run\n\n")
input('Enter')


ValQuickCommand = '''  python PrimedInclusivity -i ValQuick \
                        -v {0}/Validation/ValQuick/Target.vcf.gz \
                        -t {0}/Validation/ValQuick/Covid_19.fasta \
                        -p {0}/Validation/ValQuick/Primer_input.csv \
                        -s {0}/Validation/ValQuick/Sets_input.csv \
                        -o {0}/Validation/ValQuick/Output/ \
                        -OG {0}/pkgs/Output/Default_Output_Guide.json \
                        -RG {0}/Validation/ValQuick/RepGuide.json \
                        -OT A
                        '''.format(PackagePath)


subprocess.run(ValQuickCommand,shell=True)


RefDF = pd.read_csv('{}/Validation/ValQuick/Expected_Output/ValQuick_Sets_Table.csv'.format(PackagePath))
CheckDF = pd.read_csv('{}/Validation/ValQuick/Output/ValQuick_Sets_Table.csv'.format(PackagePath))


if assert_frame_equal(RefDF, CheckDF,check_dtype=False) == None:
    print('\n\n\nValidation Pass\n\n\n')
else:
    sys.exit()
