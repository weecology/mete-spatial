import numpy as np
import csv
import sys
import os
from math import exp

import mete

if len(sys.argv) > 1:
    S0 = int(sys.argv[1])
    N0 = int(sys.argv[2])

if os.path.exists('../demo') is False:
    os.mkdir('../demo')

beta = mete.get_beta(S0, N0)

n0 = mete.trunc_logser_rvs(exp(-beta), N0, S0)
n0 = list(n0)
n0 = [int(x) for x in n0]
n0.sort(reverse=True)

rad = mete.get_mete_rad(S0, N0)[0]

Amax = 4
Amin = 1

recur = mete.downscale_sar(Amax, S0, N0, Amin)
recur_obsSAD = mete.downscale_sar_fixed_abu(Amax, n0, Amin)

Avals = recur_obsSAD[0][ : ]

nonrecur = mete.sar_noniterative(Avals, Amax, S0, N0, 'precise')
nonrecur_obsSAD = mete.sar_noniterative_fixed_abu(Avals, Amax, n0)

sad_out = np.empty((S0, 2)) 

sad_out[ : , 0] = n0
sad_out[ : , 1] = rad

filename = '../demo/' + 'abu_sar_demo.txt'

writer = open(filename, 'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['n0', 'sad'])
for i in range(0, np.shape(sad_out)[0]):
    datawriter.writerow(sad_out[i, ])

writer.close()

sar_out = np.empty((3, 4))

sar_out[ : , 0] = recur[1] + [S0]
sar_out[ : , 1] = recur_obsSAD[1]
sar_out[ : , 2] = nonrecur[1]
sar_out[ : , 3] = nonrecur_obsSAD[1]

filename = '../demo/' + 'rich_sar_demo.txt'

writer = open(filename, 'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['recur', 'recur_obsSAD', 'nonrecur', 'nonrecur_obsSAD'])
for i in range(0, np.shape(sar_out)[0]):
    datawriter.writerow(sar_out[i, ])

writer.close()