from working_functions import *

print 'Plotting obs-pred plot for SAD, ...'

datafile = open('../data/shrtnames.txt', 'r')
datareader = csv.reader(datafile)
shrtnames = []
for row in datareader:
    shrtnames.append(row)

shrtnames = shrtnames[0][0].split()

plot_obs_pred_sad(shrtnames, data_dir='../sad/',
                  dest_dir='../figs/sad_figs/')

for shrt_name in shrtnames:
    plot_obs_pred_sad([shrt_name], data_dir='../sad/', 
                      dest_dir='../figs/sad_figs/' + shrt_name + '_')
    
print 'Plotting obs-pred plot for SAD, complete!'
