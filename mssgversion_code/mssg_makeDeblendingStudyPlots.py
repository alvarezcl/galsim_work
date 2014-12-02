# Start: 10-24-2014
# MSSG

# Takes the files on disk and plots them up 

import matplotlib.pyplot as plt
import sys
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--subdir", default="deblendsOutput/", help="output text filename")
args = parser.parse_args()

############ Read in file
fname = args.subdir + 'offset_deblendingTests_50runs.txt'

#fname = args.subdir + 'deblendingTests_50runs.txt'

fitdat = np.loadtxt(fname)

fnumvec = fitdat[:,0]

# Initze all a vecs
e1a_in = fitdat[:,1]
e1a_unbl = fitdat[:,3] #;      e1a_unblerr = []
e1a_debl = fitdat[:,4] # ;       # e1a_deblerr = []
e1a_unblresid =  e1a_unbl - e1a_in 
e1a_deblresid =  e1a_debl - e1a_in 

e1b_in = fitdat[:,5]
e1b_unbl = fitdat[:,7] #;      e1b_unblerr = []
e1b_debl = fitdat[:,8] # ;       # e1b_deblerr = []
e1b_unblresid =  e1b_unbl - e1b_in 
e1b_deblresid =  e1b_debl - e1b_in 

fnum =  fitdat[:,0]

print 'fnum = ', fnum 
print 'e1a_in = ' ,e1a_in 
print 'e1a_unbl  = ' ,e1a_unbl 
print 'e1a_debl = ', e1a_debl
print 'e1a_unblresid = ',e1a_unblresid

print 'e1b_in = ' ,e1b_in 
print 'e1b_unbl  = ' ,e1b_unbl 
print 'e1b_debl = ', e1b_debl
print 'e1b_unblresid = ',e1b_unblresid

numfiles = 50

#e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

e1a_range = [0.5,  0,  -0.5]
e1b_range = [0.5,  0,  -0.5]

xshift = [0.01, 0.01, 0.01]

e1shifted = np.array(e1a_range) + np.array(xshift)

############################################################ e1 a plots
for e1bin in e1b_range:
    e1bstr = str(e1bin)

    vece1a_in = []  
    vece1a_unbl = []
    vece1a_debl = []
    vece1a_unblresid = []  
    vece1a_deblresid = []  
    vece1a_unblerr = []  
    vece1a_deblerr = []  

    for e1ain in e1a_range:
        e1astr = str(e1ain)
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
        print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 


        vece1a_in.append( e1a_in[i].mean() )
        vece1a_unbl.append( e1a_unbl[i].mean() )
        vece1a_debl.append( e1a_debl[i].mean() )
        vece1a_unblresid.append (e1a_unblresid[i].mean()  )
        vece1a_deblresid.append (e1a_deblresid[i].mean() )
        vece1a_unblerr.append( e1a_unbl[i].std()  )
        vece1a_deblerr.append( e1a_debl[i].std() ) 

######## e1 a plots
    print 'vece1a_unbl = ', vece1a_unbl
    plt.figure(figsize=(15,12))
    xlimit = 0.6;    ylimit = 0.15
    plt.xlim( -xlimit, xlimit);    plt.ylim( -ylimit, ylimit)
   
    plt.scatter( e1shifted , vece1a_unblresid, color = 'g' , s=50.0  )
    gline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1a" , linewidth=4.0)
    yline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr/np.sqrt(numfiles),  ecolor='y',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

    plt.scatter( e1a_range, vece1a_deblresid, color = 'b' , s=50.0  )
    bline = plt.errorbar(e1a_range, vece1a_deblresid, vece1a_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1a" , linewidth= 4.0)
    rline = plt.errorbar(e1a_range, vece1a_deblresid, vece1a_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit, error/sqrt(N)" , linewidth= 8.0 )

    plt.title("Resids for fits with e1bin = " +e1bstr )
    plt.legend() # (handles=[gline,bline])
    plt.xlabel('$e_{1a in}$',fontsize=18)
    plt.ylabel('$e_{1a fit}-e_{1a in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
    plt.show()

#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")

#e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

e1a_range = [0.5, 0, -0.5]
e1b_range = [0.5, 0, -0.5]



############################################################ e1 b plots
for e1ain in e1a_range:
    e1astr = str(e1ain)

    vece1b_in = []  
    vece1b_unbl = []
    vece1b_debl = []
    vece1b_unblresid = []  
    vece1b_deblresid = []  
    vece1b_unblerr = []  
    vece1b_deblerr = []  

    for e1bin in e1b_range:
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
        print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        vece1b_in.append( e1b_in[i].mean() )
        vece1b_unbl.append( e1b_unbl[i].mean() )
        vece1b_debl.append( e1b_debl[i].mean() )
        vece1b_unblresid.append (e1b_unblresid[i].mean()  )
        vece1b_deblresid.append (e1b_deblresid[i].mean() )
        vece1b_unblerr.append( e1b_unbl[i].std()  )
        vece1b_deblerr.append( e1b_debl[i].std() ) 

######## e1 b plots
    print 'vece1b_unbl = ', vece1b_unbl
    plt.figure(figsize=(15,12))
    xlimit = 0.6;    ylimit = 0.15
    plt.xlim( -xlimit, xlimit);    plt.ylim( -ylimit, ylimit)
   
    plt.scatter( e1shifted , vece1b_unblresid, color = 'g' , s=50.0  )
    gline = plt.errorbar( e1shifted, vece1b_unblresid, vece1b_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1b" , linewidth=4.0)
    yline = plt.errorbar( e1shifted, vece1b_unblresid, vece1b_unblerr/np.sqrt(numfiles),  ecolor='y',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

    plt.scatter( e1b_range, vece1b_deblresid, color = 'b' , s=50.0  )
    bline = plt.errorbar(e1b_range, vece1b_deblresid, vece1b_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1b" , linewidth= 4.0)
    rline = plt.errorbar(e1b_range, vece1b_deblresid, vece1b_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit, error/sqrt(N)" , linewidth= 8.0 )

    plt.title("Resids for fits with e1ain = " +e1astr )
    plt.legend()
    plt.xlabel('$e_{1b in}$',fontsize=18)
    plt.ylabel('$e_{1b fit}-e_{1b in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
    plt.show()

sys.exit()
