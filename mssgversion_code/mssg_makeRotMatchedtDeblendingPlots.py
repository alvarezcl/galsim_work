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
fname = args.subdir + 'deblendingTests_50runs.txt' # Orig


# fname = args.subdir + 'offsetEachQuarterPixelAwayFromCenter_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBothQuarterPixelLeft_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBoth0.005PixelLeft_deblendingTests_50runs.txt'     # 
#fname = args.subdir + 'offsetBoth0.005PixelRight_deblendingTests_50runs.txt'


#################################### Load up orig data
#fname = args.subdir + 'deblendingTests_peak_A_(-1, 0)__peak_B_(1, 0)_5_runs.txt' # Orig
fname = args.subdir + 'deblendingTests_50runs.txt' # Orig
   # Horiz sep- move centers to -eps and +eps
#fname = args.subdir + 'offsetLeftOne200thPixelLeftRightOne200thPixelRight_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to +eps and -eps
# fname = args.subdir + 'offsetLeftOne200thPixelRightRightOne200thPixelLeft_deblendingTests_50runs.txt' 

# Load file
fitdat = np.loadtxt(fname)

# File num
fnum =  fitdat[:,0]

try:
# Initze all a vecs
    orige1a_in = fitdat[:,1]
    orige2a_in = fitdat[:,2]
    orige1a_unbl = fitdat[:,3] 
    orige1a_debl = fitdat[:,4] 
    orige2a_unbl = fitdat[:,5] 
    orige2a_debl = fitdat[:,6] 
    orige1a_unblresid =  orige1a_unbl - orige1a_in 
    orige1a_deblresid =  orige1a_debl - orige1a_in 
    orige2a_unblresid =  orige2a_unbl - orige2a_in 
    orige2a_deblresid =  orige2a_debl - orige2a_in 
    
    orige1b_in = fitdat[:,7]
    orige2b_in = fitdat[:,8]
    orige1b_unbl = fitdat[:,9] 
    orige1b_debl = fitdat[:,10]
    orige2b_unbl = fitdat[:,11] 
    orige2b_debl = fitdat[:,12]
    orige1b_unblresid =  orige1b_unbl - orige1b_in 
    orige1b_deblresid =  orige1b_debl - orige1b_in 
    orige2b_unblresid =  orige2b_unbl - orige2b_in 
    orige2b_deblresid =  orige2b_debl - orige2b_in 

except:

    print "*********************** No e2 in output file, reading just the orig files"
    orige1a_in = fitdat[:,1]
    orige2a_in = fitdat[:,2]
    orige1a_unbl = fitdat[:,3] 
    orige1a_debl = fitdat[:,4] 
    orige1a_unblresid =  orige1a_unbl - orige1a_in 
    orige1a_deblresid =  orige1a_debl - orige1a_in 

    orige1b_in = fitdat[:,5]
    orige2b_in = fitdat[:,6]
    orige1b_unbl = fitdat[:,7] 
    orige1b_debl = fitdat[:,8]
    orige1b_unblresid =  orige1b_unbl - orige1b_in 
    orige1b_deblresid =  orige1b_debl - orige1b_in 


#################################### Load up nonrot data
# Horiz sep- Random half pixel offset for both
fname = args.subdir + 'deblendingTests_peak_A_(-1, 0)__peak_B_(1, 0)_50_runsAndRandomOffsetHalfPixelEach.txt'

fitdat = np.loadtxt(fname)

# File num
fnum =  fitdat[:,0]

# Initze all a vecs
e1a_in = fitdat[:,1]
e2a_in = fitdat[:,2]
e1a_unbl = fitdat[:,3] 
e1a_debl = fitdat[:,4] 
e2a_unbl = fitdat[:,5] 
e2a_debl = fitdat[:,6] 
e1a_unblresid =  e1a_unbl - e1a_in 
e1a_deblresid =  e1a_debl - e1a_in 
e2a_unblresid =  e2a_unbl - e2a_in 
e2a_deblresid =  e2a_debl - e2a_in 

e1b_in = fitdat[:,7]
e2b_in = fitdat[:,8]
e1b_unbl = fitdat[:,9] 
e1b_debl = fitdat[:,10]
e2b_unbl = fitdat[:,11] 
e2b_debl = fitdat[:,12]
e1b_unblresid =  e1b_unbl - e1b_in 
e1b_deblresid =  e1b_debl - e1b_in 
e2b_unblresid =  e2b_unbl - e2b_in 
e2b_deblresid =  e2b_debl - e2b_in 

#################################### Load up rotated data 
rotfname = args.subdir + 'deblendingTests_peak_A_(0, -1)__peak_B_(0, 1)_50_runsAndRandomOffsetHalfPixelEach.txt'

fitdat = np.loadtxt(rotfname)

# File num
rotfnum =  fitdat[:,0]

# Initze all a vecs
rote1a_in = fitdat[:,1]
rote2a_in = fitdat[:,2]
rote1a_unbl = fitdat[:,3] 
rote1a_debl = fitdat[:,4] 
rote2a_unbl = fitdat[:,5] 
rote2a_debl = fitdat[:,6] 
rote1a_unblresid =  rote1a_unbl - rote1a_in 
rote1a_deblresid =  rote1a_debl - rote1a_in 
rote2a_unblresid =  rote2a_unbl - rote2a_in 
rote2a_deblresid =  rote2a_debl - rote2a_in 

rote1b_in = fitdat[:,7]
rote2b_in = fitdat[:,8]
rote1b_unbl = fitdat[:,9] 
rote1b_debl = fitdat[:,10]
rote2b_unbl = fitdat[:,11] 
rote2b_debl = fitdat[:,12]
rote1b_unblresid =  rote1b_unbl - rote1b_in 
rote1b_deblresid =  rote1b_debl - rote1b_in 
rote2b_unblresid =  rote2b_unbl - rote2b_in 
rote2b_deblresid =  rote2b_debl - rote2b_in 

# For printing
'''
print 'rotfnum = ', rotfnum 
print 'rote1a_in = ' ,rote1a_in 
print 'rote1a_unbl  = ' ,rote1a_unbl 
print 'rote1a_debl = ', rote1a_debl
print 'rote1a_unblresid = ',rote1a_unblresid

print 'rote1b_in = ' ,rote1b_in 
print 'rote1b_unbl  = ' ,rote1b_unbl 
print 'rote1b_debl = ', rote1b_debl
print 'rote1b_unblresid = ',rote1b_unblresid
'''

################### Now declare some needed  stuff
numfiles = 50

## We know the inputs were made at these intervals
e1a_range = [0.5,  0,  -0.5]
e1b_range = [0.5,  0,  -0.5]

## This is when we want to plot several vecs at the same input xval, we put in a shift by hand
xshift = [0.01, 0.01, 0.01]
xshiftL = [-0.01, -0.01, -0.01]
xshiftLL = [-0.02, -0.02, -0.02]
xshiftRR = [+0.02, 0.02, 0.02]

e1shifted = np.array(e1a_range) + np.array(xshift)
e1Lshifted = np.array(e1a_range) + np.array(xshiftL)
e1LLshifted = np.array(e1a_range) + np.array(xshiftLL)
e1RRshifted = np.array(e1a_range) + np.array(xshiftRR)

## If we're drawing histos we need this
nbins = 10 # Set num of bins for histo

############################################################ e1 a plots
for e1bin in e1b_range:
    e1bstr = str(e1bin)

    # Declare vecs for nonrot file
    vece1a_in = []  
    vece1a_unbl = []
    vece1a_debl = []
    vece1a_unblresid = []  
    vece1a_deblresid = []  
    vece1a_unblerr = []  
    vece1a_deblerr = []  

    vece2a_in = []  
    vece2a_unbl = []
    vece2a_debl = []
    vece2a_unblresid = []  
    vece2a_deblresid = []  
    vece2a_unblerr = []  
    vece2a_deblerr = []  

    # Declare vecs for rot file
    rotvece1a_in = []  
    rotvece1a_unbl = []
    rotvece1a_debl = []
    rotvece1a_unblresid = []  
    rotvece1a_deblresid = []  
    rotvece1a_unblerr = []  
    rotvece1a_deblerr = []  

    rotvece2a_in = []  
    rotvece2a_unbl = []
    rotvece2a_debl = []
    rotvece2a_unblresid = []  
    rotvece2a_deblresid = []  
    rotvece2a_unblerr = []  
    rotvece2a_deblerr = []  

    # Declare vecs for orig file
    origvece1a_in = []  
    origvece1a_unbl = []
    origvece1a_debl = []
    origvece1a_unblresid = []  
    origvece1a_deblresid = []  
    origvece1a_unblerr = []  
    origvece1a_deblerr = []  

    origvece2a_in = []  
    origvece2a_unbl = []
    origvece2a_debl = []
    origvece2a_unblresid = []  
    origvece2a_deblresid = []  
    origvece2a_unblerr = []  
    origvece2a_deblerr = []  

    #################### Run over all e1a fits
    for e1ain in e1a_range:
        e1astr = str(e1ain)
        # Get indices for all runs that match the inputs for the nonrot file
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
#        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
 #       print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        # Get indices for all runs that match the inputs for the nonrot file for the rotated file
        j = np.where(np.logical_and(rote1a_in == -e1ain, rote1b_in == -e1bin))  # When we rotate by 90 deg the sign of e1 flips

        # Get indices for all runs that match the inputs for the nonrot file for the orig file
        k = np.where(np.logical_and(orige1a_in == e1ain, orige1b_in == e1bin))  

        # Take avg and stddev for all vecs from nonrot file (these vecs will be just number of e1ain vals long -- i.e. 3 for [-0.5, 0, 0.5])
        vece1a_in.append( e1a_in[i].mean() )
        vece1a_unbl.append( e1a_unbl[i].mean() )
        vece1a_debl.append( e1a_debl[i].mean() )
        vece1a_unblresid.append (e1a_unblresid[i].mean()  )
        vece1a_deblresid.append (e1a_deblresid[i].mean() )
        vece1a_unblerr.append( e1a_unbl[i].std()  )
        vece1a_deblerr.append( e1a_debl[i].std() ) 

        vece2a_in.append( e2a_in[i].mean() )
        vece2a_unbl.append( e2a_unbl[i].mean() )
        vece2a_debl.append( e2a_debl[i].mean() )
        vece2a_unblresid.append (e2a_unblresid[i].mean()  )
        vece2a_deblresid.append (e2a_deblresid[i].mean() )
        vece2a_unblerr.append( e2a_unbl[i].std()  )
        vece2a_deblerr.append( e2a_debl[i].std() ) 

        # Take avg and stddev for all vecs from rot file
        rotvece1a_in.append( rote1a_in[j].mean() )
        rotvece1a_unbl.append( rote1a_unbl[j].mean() )
        rotvece1a_debl.append( rote1a_debl[j].mean() )
        rotvece1a_unblresid.append ( rote1a_unblresid[j].mean()  )
        rotvece1a_deblresid.append ( rote1a_deblresid[j].mean() )
        rotvece1a_unblerr.append( rote1a_unbl[j].std()  )
        rotvece1a_deblerr.append( rote1a_debl[j].std() ) 

        rotvece2a_in.append( rote2a_in[j].mean() )
        rotvece2a_unbl.append( rote2a_unbl[j].mean() )
        rotvece2a_debl.append( rote2a_debl[j].mean() )
        rotvece2a_unblresid.append ( rote2a_unblresid[j].mean()  )
        rotvece2a_deblresid.append ( rote2a_deblresid[j].mean() )
        rotvece2a_unblerr.append( rote2a_unbl[j].std()  )
        rotvece2a_deblerr.append( rote2a_debl[j].std() ) 

        # Take avg and stddev for all vecs from orig file (these vecs will be just number of e1ain vals long -- i.e. 3 for [-0.5, 0, 0.5])
        origvece1a_in.append( orige1a_in[k].mean() )
        origvece1a_unbl.append( orige1a_unbl[k].mean() )
        origvece1a_debl.append( orige1a_debl[k].mean() )
        origvece1a_unblresid.append (e1a_unblresid[k].mean()  )
        origvece1a_deblresid.append (e1a_deblresid[k].mean() )
        origvece1a_unblerr.append( orige1a_unbl[k].std()  )
        origvece1a_deblerr.append( orige1a_debl[k].std() ) 

        origvece2a_in.append( orige2a_in[k].mean() )
        origvece2a_unbl.append( orige2a_unbl[k].mean() )
        origvece2a_debl.append( orige2a_debl[k].mean() )
        origvece2a_unblresid.append (e2a_unblresid[k].mean()  )
        origvece2a_deblresid.append (e2a_deblresid[k].mean() )
        origvece2a_unblerr.append( orige2a_unbl[k].std()  )
        origvece2a_deblerr.append( orige2a_debl[k].std() ) 

        '''
        #  Make histo of e1a fit vals
        plt.title("Histo of e1a debl fit dist for e1ain = " +e1astr + " with e1bin = " +e1bstr )
        plt.hist(e1a_debl[i],bins= nbins)
        plt.show()
        '''

#    sys.exit()

######## e1 a plots
    print 'vece1a_unbl = ', vece1a_unbl

    print " Using file ", fname
    plt.figure(figsize=(15,12))
    xlimit = 0.6;    ylimit = 0.15
#    plt.xlim( -xlimit, xlimit);    plt.ylim( -ylimit, ylimit)

#### Unblended fit e1a plots
# (Note we are horizontally offsetting these points by xshift, as defined above)
    plt.scatter( e1shifted , vece1a_unblresid, color = 'g' , s=50.0  ) # The central value point -- 's' here is the size of the point -- 
    gline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1a" , linewidth=4.0)
    yline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit e1a plots
    plt.scatter( e1a_range, vece1a_deblresid, color = 'b' , s=50.0  )
    bline = plt.errorbar(e1a_range, vece1a_deblresid, vece1a_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1a" , linewidth= 4.0)
    rline = plt.errorbar(e1a_range, vece1a_deblresid, vece1a_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit rote1a plots
    plt.scatter( e1Lshifted, rotvece1a_deblresid, color = 'm' , s=50.0  )
    bline = plt.errorbar(e1Lshifted, rotvece1a_deblresid, rotvece1a_deblerr,  ecolor='c',linestyle=' ', label = "Resid for deblended fit for rotvece1a" , linewidth= 4.0)
    rline = plt.errorbar(e1Lshifted, rotvece1a_deblresid, rotvece1a_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit for rotvece1a, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit orig plots
    plt.scatter( e1RRshifted, origvece1a_deblresid, color = 'y' , s=50.0  )
    oline = plt.errorbar(e1RRshifted, origvece1a_deblresid, origvece1a_deblerr,  ecolor='y',linestyle=' ', label = "Resid for deblended fit for origvece1a" , linewidth= 4.0)
    ooline = plt.errorbar(e1RRshifted, origvece1a_deblresid, origvece1a_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit for origvece1a, error/sqrt(N)" , linewidth= 8.0 )


#### Deblended fit e1a - rote1a plots
 
    e1sumvec = np.array(vece1a_deblresid) + np.array(rotvece1a_deblresid) 
    e1errvec =  np.sqrt(np.array(vece1a_deblerr)**2 + np.array(rotvece1a_deblerr)**2)

    plt.scatter( e1LLshifted, e1sumvec , color = 'c' , s=50.0  )
    bline = plt.errorbar(e1LLshifted, e1sumvec , e1errvec,  ecolor='black',linestyle=' ', label = "Deblended fit - rotated fit, for e1a  " , linewidth= 4.0)
    rline = plt.errorbar(e1LLshifted, e1sumvec , e1errvec/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Deblended fit - rotated fit, for e1a, error/sqrt(N)  " , linewidth= 8.0)


    print "vece1a_deblresid =  ", vece1a_deblresid 
    print "rotvece1a_deblresid =  ", rotvece1a_deblresid 



    '''
#### Unblended e2a fit plots
    plt.scatter( e1LLshifted, vece2a_unblresid, color = 'c' , s=50.0  )
    bline = plt.errorbar(e1LLshifted, vece2a_unblresid, vece2a_unblerr,  ecolor='black',linestyle=' ', label = "Resid for unblended fit for e2a" , linewidth= 4.0)
    rline = plt.errorbar(e1LLshifted, vece2a_unblresid, vece2a_unblerr/np.sqrt(numfiles),  ecolor='m',linestyle=' ', label = "Resid for unblended e2a fit, error/sqrt(N)" , linewidth= 8.0 )


#### Deblended e2a fit plots
    plt.scatter( e1Lshifted, vece2a_deblresid, color = 'c' , s=50.0  )
    bline = plt.errorbar(e1Lshifted, vece2a_deblresid, vece2a_deblerr,  ecolor='c',linestyle=' ', label = "Resid for deblended fit for e2a" , linewidth= 4.0)
    rline = plt.errorbar(e1Lshifted, vece2a_deblresid, vece2a_deblerr/np.sqrt(numfiles),  ecolor='m',linestyle=' ', label = "Resid for deblended e2a fit, error/sqrt(N)" , linewidth= 8.0 )
    '''
## Plot title and legend
    plt.title("Resids for fits with e1bin = " +e1bstr )
    plt.legend() # (handles=[gline,bline])
    plt.xlabel('$e_{1a in}$',fontsize=18)
    plt.ylabel('$e_{1a fit}-e_{1a in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
    plt.show()

#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")

############################################################ e1 b plots
for e1ain in e1a_range:
    e1astr = str(e1ain)

    # Declare vecs for nonrot file
    vece1b_in = []  
    vece1b_unbl = []
    vece1b_debl = []
    vece1b_unblresid = []  
    vece1b_deblresid = []  
    vece1b_unblerr = []  
    vece1b_deblerr = []  

    # Declare vecs for rot file
    rotvece1b_in = []  
    rotvece1b_unbl = []
    rotvece1b_debl = []
    rotvece1b_unblresid = []  
    rotvece1b_deblresid = []  
    rotvece1b_unblerr = []  
    rotvece1b_deblerr = []  

    for e1bin in e1b_range:
        # Get indices for all runs that match the inputs for the nonrot file
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
#        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
 #       print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        # Get indices for all runs that match the inputs for the nonrot file for the rotated file
        j = np.where(np.logical_and(rote1a_in == -e1ain, rote1b_in == -e1bin))  # When we rotate by 90 deg the sign of e1 flips

        # Take avg and stddev for all vecs from nonrot file (these vecs will be just number of e1ain vals long -- i.e. 3 for [-0.5, 0, 0.5])
        vece1b_in.append( e1b_in[i].mean() )
        vece1b_unbl.append( e1b_unbl[i].mean() )
        vece1b_debl.append( e1b_debl[i].mean() )
        vece1b_unblresid.append (e1b_unblresid[i].mean()  )
        vece1b_deblresid.append (e1b_deblresid[i].mean() )
        vece1b_unblerr.append( e1b_unbl[i].std()  )
        vece1b_deblerr.append( e1b_debl[i].std() ) 

        # Take avg and stddev for all vecs from rot file
        rotvece1b_in.append( rote1b_in[j].mean() )
        rotvece1b_unbl.append( rote1b_unbl[j].mean() )
        rotvece1b_debl.append( rote1b_debl[j].mean() )
        rotvece1b_unblresid.append ( rote1b_unblresid[j].mean()  )
        rotvece1b_deblresid.append ( rote1b_deblresid[j].mean() )
        rotvece1b_unblerr.append( rote1b_unbl[j].std()  )
        rotvece1b_deblerr.append( rote1b_debl[j].std() ) 


        '''
        #  Make histo of e1a fit vals
        plt.title("Histo of e1b debl fit dist for e1ain = " +e1astr + " with e1bin = " +e1bstr )
        plt.hist(e1b_debl[i],bins= nbins)
        plt.show()
        '''

######## e1 b plots
    print 'vece1b_unbl = ', vece1b_unbl
    plt.figure(figsize=(15,12))
    xlimit = 0.6;    ylimit = 0.15
#    plt.xlim( -xlimit, xlimit);    plt.ylim( -ylimit, ylimit)

#### Unblended fit e1b plots
# (Note we are horizontally offsetting these points by xshift, as defined above)   
    plt.scatter( e1shifted , vece1b_unblresid, color = 'g' , s=50.0  )
    gline = plt.errorbar( e1shifted, vece1b_unblresid, vece1b_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1b" , linewidth=4.0)
    yline = plt.errorbar( e1shifted, vece1b_unblresid, vece1b_unblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit e1b plots
    plt.scatter( e1b_range, vece1b_deblresid, color = 'b' , s=50.0  )
    bline = plt.errorbar(e1b_range, vece1b_deblresid, vece1b_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1b" , linewidth= 4.0)
    rline = plt.errorbar(e1b_range, vece1b_deblresid, vece1b_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit rote1b plots
    plt.scatter( e1Lshifted, rotvece1b_deblresid, color = 'm' , s=50.0  )
    bline = plt.errorbar(e1Lshifted, rotvece1b_deblresid, rotvece1b_deblerr,  ecolor='c',linestyle=' ', label = "Resid for deblended fit for rotvece1b" , linewidth= 4.0)
    rline = plt.errorbar(e1Lshifted, rotvece1b_deblresid, rotvece1b_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit for rotvece1b, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit e1b - rote1b plots
    e1sumvec = np.array(vece1b_deblresid) + np.array(rotvece1b_deblresid) 
    e1errvec =  np.sqrt(np.array(vece1b_deblerr)**2 + np.array(rotvece1b_deblerr)**2)

    plt.scatter( e1LLshifted, e1sumvec , color = 'c' , s=50.0  )
    bline = plt.errorbar(e1LLshifted, e1sumvec , e1errvec,  ecolor='black',linestyle=' ', label = "Deblended fit - rotated fit, for e1b  " , linewidth= 4.0)
    rline = plt.errorbar(e1LLshifted, e1sumvec , e1errvec/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Deblended fit - rotated fit, for e1b, error/sqrt(N)  " , linewidth= 8.0)

## Plot title and legend
    plt.title("Resids for fits with e1ain = " +e1astr )
    plt.legend()
    plt.xlabel('$e_{1b in}$',fontsize=18)
    plt.ylabel('$e_{1b fit}-e_{1b in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
    plt.show()

sys.exit()


# Other colors to use: y = yellow, m = magenta
