"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the universal curve
"""

import csv
import sys
import os

import mete


data_names = ['bci','cocoli1','cocoli2','cross','serp','sherman1','sherman2','sherman3']
slopes = []
for shrt_name in data_names:
    datafile = open('./sar/' + shrt_name + '_empir_sar.csv', 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
   
    site_data = []
    for i in range(1,len(data)):
        site_data.append(data[i][0:3])

    site_data = [map(float,x) for x in site_data]

    slopes.append(mete.get_slopes(site_data))

writer = open('./sar/sar_slopes.csv','wb') 
datawriter = csv.writer(writer)
datawriter.writerow(['comm','area', 'obs z', 'pred z', 'N/S'])
for i in range(0,len(data_names)):
    for j in range(0,len(slopes[i])):
        slopes[i][j].insert(0,data_names[i])
        datawriter.writerow(slopes[i][j])

writer.close()  


