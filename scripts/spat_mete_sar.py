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

site_names = ['bci','cocoli1','cocoli2','cross','sherman1','sherman2',
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak']

for shrt_name in site_names:
    datafile = open('../sar/' + shrt_name + '_empir_sar.csv', 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
   
    # drop the header row
    site_data = []
    for i in range(1,len(data)):
        site_data.append(data[i][0:3])

    site_data = [map(float,x) for x in site_data]

    # enforce a minimum individual density of 2
    indices = mete.which([site_data[i][2] > 2 for i in range(0, len(site_data))])

    site_data = [site_data[i] for i in indices]

    site_data = np.array(site_data)

    # get parameters needed for computing the mete sar
    Amin = min(site_data[:, 0])
    Amax = max(site_data[:, 0])
    S = max(site_data[:, 1])
    N = max(site_data[:, 2])

    sar_down = mete.downscale_sar(Amax, S, N, Amin)

    # add values at Amax
    sar_down[0].append(Amax)
    sar_down[1].append(S)

    # Make an array so that the data is easier to output
    out = np.empty((len(sar_down[0]), 2)) 
    for i in range(0,2):
        out[ : , i] = sar_down[i] 

    filename = '../sar/' + shrt_name + '_mete_sar.txt'
    writer = open(filename, 'wb') 
    datawriter = csv.writer(writer)
    datawriter.writerow(['area', 'sr'])
    for i in range(0, len(sar_down[0])):
        datawriter.writerow(out[i, ])

    writer.close()

