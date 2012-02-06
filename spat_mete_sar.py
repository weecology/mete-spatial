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

if 'sar' not in os.listdir(os.path.curdir):
    os.mkdir('sar')

if(len(sys.argv) > 1):
    S = int(sys.argv[1]) 
    N = int(sys.argv[2]) 
    bisec = int(sys.argv[4])
    shrt_name = sys.argv[5]
else:
    S = 10
    N = 100
    bisec = 4
    shrt_name = 'test'
    
Amax = 2 ** (bisec - 1) # number of quadrats per community 
Amin = 1
 
sar_down = mete.downscale_sar(Amax,S,N,Amin)

# Make an array so that the data is easier to output
out = np.empty((bisec-1, 2)) 
for i in range(0,2):
    out[:,i] = sar_down[i]

filename = './sar/' + shrt_name + '_mete_sar.txt'

writer = open(filename,'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['area','sr'])
for i in range(0, bisec-1):
    datawriter.writerow(out[i,])

writer.close()  

