from General.Core import *
from shutil import copyfile
import filecmp

'''
#https://primer3.ut.ee/primer3web_help.htm#PRIMER_TM_FORMULA
Primer3Primer = 'CGTGACGTGACGGACT'

# Delta H  == -128.8 kcal/mol
# Delta S  == -345.2 cal/k*mol

x = Thermodynamics(Primer3Primer,'GCACTGCACTGCCTGA')
print(x)
'''


################# Manual Example
ME = 'ATTAAAGGTTTATACCTTCCCAGGTAAC'
y = Thermodynamics('ATTAAAGGTTTATACCTTCCCAGGTAAC',
                   'TAATTTCCAAATATAGAAGGGTCCATTG')

print(y)


ReagentConcs = {'DivalentSaltConc':2,
                 'MonovalentSaltConc':50,
                 'dNTPConc':0.8,
                 'OligoConc':0.5,
                 'RunningTm':65}


print(ThermoSaltCorrection(y,
				     (len(ME)*2)-2,
					 ReagentConcs,
					 RunningTm=65))
