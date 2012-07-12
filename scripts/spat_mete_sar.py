"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the SAR
"""

from math import log
import numpy as np
import csv
import sys
import os

import mete

if 'sar' not in os.listdir('..'):
    os.mkdir('../sar')

shrtname = sys.argv[1]
    
datafile = open('../sar/empir_sars.csv', 'r')
datareader = csv.reader(datafile)
data = []
for row in datareader:
    data.append(row)

for i in range(1, len(data)):
    data[i][1 : ] = map(float, data[i][1 : ])

indices = mete.which([data[i][0] == shrtname and data[i][3] > 2 for i in range(0, len(data))])

Amin = min([int(data[i][1]) for i in indices])
Amax = max([int(data[i][1]) for i in indices])
S = max([int(data[i][2]) for i in indices])
N = max([int(data[i][3]) for i in indices])

sar_down = mete.downscale_sar(Amax, S, N, Amin)

# add values at Amax
sar_down[0].append(Amax)
sar_down[1].append(S)

# Make an array so that the data is easier to output

out = np.empty((len(sar_down[0]), 2)) 
for i in range(0,2):
    out[ : , i] = sar_down[i] 

filename = '../sar/' + shrtname + '_mete_sar.txt'

writer = open(filename, 'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['area', 'sr'])
for i in range(0, len(sar_down[0])):
    datawriter.writerow(out[i, ])

writer.close()  

