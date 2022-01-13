import os
import json
import subprocess
import pandas as pd
import filecmp
from pandas._testing import assert_frame_equal

# Extract location of the PrimedInclusivity Software
PackagePath = os.path.dirname(os.path.abspath(__file__))


print("Tutorial Example")
input('Enter:')

PIDefaultCustomCommand = 'python PIConfig -D'
subprocess.run(PIDefaultCustomCommand,shell=True)


# Create Custom RPA Config Json File
CustomJson = {"PIConfig":{
                         "Classification":['{}/pkgs/Classification/DefaultClassification.py'.format(PackagePath),
                                           '{}/Tutorial/MinLenClassification.py'.format(PackagePath)],

                         "Output":['{}/pkgs/Output/DefaultOutput.py'.format(PackagePath),
                                   '{}/Tutorial/MinLenOutput.py'.format(PackagePath)],
                            }}

ConfigPath = '{}/Tutorial/Tutorial_PIConfigVal.json'.format(PackagePath)

with open(ConfigPath, 'w') as CJP:
    json.dump(CustomJson, CJP)
CJP.close()



ConfigCommand = 'python PIConfig -CF {}'.format(ConfigPath)
subprocess.run(ConfigCommand,shell=True)
print('Tutorial Configuration Complete')



TutorialCommand = ''' python PrimedInclusivity -i Tutorial_Screen \
                    -v {0}/Tutorial/Target.vcf.gz \
                    -t {0}/Tutorial/VirusA.fasta \
                    -p {0}/Tutorial/MinLen_Primer_input.csv \
                    -s {0}/Tutorial/MinLen_Sets_input.csv \
                    -o {0}/Tutorial/Output/ \
                    -OG {0}/Tutorial/MinLen_Output_Guide.json \
                    -OT A
                    '''.format(PackagePath)



subprocess.run(TutorialCommand,shell=True)
