"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the universal curve
"""

import csv
import sys
import os

import mete


site_names = ['bci','cocoli1','cocoli2','cross','sherman1','sherman2',
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak']

slopes = []
for shrt_name in site_names:
    datafile = open('../sar/' + shrt_name + '_empir_sar.csv', 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
   
    site_data = []
    for i in range(1, len(data)):
        site_data.append(data[i][0:3])

    site_data = [map(float,x) for x in site_data]

    slopes.append(mete.get_slopes(site_data))

    print shrt_name

writer = open('../sar/sar_slopes.csv', 'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['comm', 'area', 'obsZ', 'predZ', 'NS'])
for i in range(0, len(site_names)):
    for j in range(0, len(slopes[i])):
        slopes[i][j].insert(0, site_names[i])
        datawriter.writerow(slopes[i][j])

writer.close()
