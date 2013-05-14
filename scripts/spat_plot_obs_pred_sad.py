from working_functions import *

print 'Plotting obs-pred plot for SAD, ...'

datasets = ['bci','cocoli1','cocoli2','cross','sherman1','sherman2',
            'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
            'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
            'bigoak']

plot_obs_pred_sad(datasets, data_dir='../sad/',
                  dest_dir='../figs/sad_figs/')

for shrt_name in datasets:
    plot_obs_pred_sad([shrt_name], data_dir='../sad/', 
                      dest_dir='../figs/sad_figs/' + shrt_name + '_')
    
print 'Plotting obs-pred plot for SAD, complete!'
