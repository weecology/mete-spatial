"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the SAR
"""

import numpy as np
import csv
import sys
import os

import mete

if 'sar' not in os.listdir('..'):
    os.mkdir('../sar')

if(len(sys.argv) > 1):
    shrt_name = sys.argv[1]
else:
    S = 10
    N = 100
    bisec = 4
    shrt_name = 'test'
    
datafile = open('../sar/empir_sars.csv', 'r')
datareader = csv.reader(datafile)
data = []
for row in datareader:
    data.append(row)

for i in range(1, len(data)):
    data[i][1 : ] = map(float, data[i][1 : ])

indices = mete.which([data[i][0] == shrtname and data[i][3] > 1 for i in range(0, len(data))])

Amin = min([int(data[i][1]) for i in indices])
Amax = max([int(data[i][1]) for i in indices])
S = max([int(data[i][2]) for i in indices])
N = max([int(data[i][3]) for i in indices])

sar_down = mete.downscale_sar(Amax, S, N, Amin)

# Make an array so that the data is easier to output
out = np.empty((bisec, 2)) 
for i in range(0,2):
    out[:,i] = sar_down[i]

filename = '../sar/' + shrt_name + '_mete_sar.txt'

writer = open(filename,'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['area','sr'])
for i in range(0, bisec):
    datawriter.writerow(out[i,])

writer.close()  

