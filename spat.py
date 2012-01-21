
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def comm_filename(S, N, ncomm, bisec, transect = False, abu = None,
                   shrt_name = None):
    """Creates the file name of a simulated community produced by the function
    sim_spatial_whole() in the mete.py module
    
    Keyword arguments:
    S -- the number of species
    N -- the number of individuals
    ncomm -- the number of communities that were generated
    bisec -- the number of bisections
    transect -- a boolean if False then it is assumed a 2D grid was generated
    abu -- the path to an abundance file else this should be None
    shrt_name -- the short hand name for the output community, this overrides
                 S and N in the output file name.
    """
                 
    if(transect):
        if(comm_name is None):  
            if(abu is None):
                filename = ('./comms/simulated_comms_S' + str(S) + '_N' +
                            str(N) + '_C' + str(ncomm) + '_B' + str(bisec) + 
                            '_transect.txt')
            else:
                filename = ('./comms/simulated_comms_S' + str(S) + '_N' +
                            str(N) + '_empirSAD' + '_C' +  str(ncomm) + '_B' +
                            str(bisec) + '_transect.txt')
        else:
            if(abu is None):
                filename = ('./comms/simulated_comms_' + comm_name + '_C' +
                            str(ncomm) + '_B' + str(bisec) + '_transect.txt')
            else:
                filename = ('./comms/simulated_comms_' + comm_name + '_empirSAD'
                            + '_C' + str(ncomm) + '_B' + str(bisec) + 
                            '_transect.txt')
    else:
        if(comm_name is None):
            if(abu is None):
                filename = ('./comms/simulated_comms_S' + str(S) + '_N' +
                            str(N) + '_C' + str(ncomm) + '_B' + str(bisec) +
                            '_grid.txt')
            else:
                filename = ('./comms/simulated_comms_S' + str(S) + '_N' + 
                            str(N) + '_empirSAD' + '_C' + str(ncomm) + '_B' +
                            str(bisec) + '_grid.txt')
        else:
            if(abu is None):
                filename = ('./comms/simulated_comms_' + comm_name + '_C' +
                            str(ncomm) + '_B' + str(bisec) + '_grid.txt')
            else:
                filename = ('./comms/simulated_comms_' + comm_name + '_empirSAD'
                            + '_C' + str(ncomm) + '_B' + str(bisec) + '_grid.txt')
    return(filename)
