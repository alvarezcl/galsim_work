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
# fname = args.subdir + 'deblendingTests_50runs.txt' # Orig
# fname = args.subdir + 'offsetEachQuarterPixelAwayFromCenter_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBothQuarterPixelLeft_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBoth0.005PixelLeft_deblendingTests_50runs.txt'     # 
#fname = args.subdir + 'offsetBoth0.005PixelRight_deblendingTests_50runs.txt'


   # Horiz sep- move centers to -eps and -eps
# fname = args.subdir + 'offsetBothOne200thPixelLeft_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to -eps and +eps
#fname = args.subdir + 'offsetLeftOne200thPixelLeftRightOne200thPixelRight_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to +eps and -eps
#fname = args.subdir + 'offsetLeftOne200thPixelRightRightOne200thPixelLeft_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to +eps and +eps
#fname = args.subdir + 'offsetBothOne200thPixelRight_deblendingTests_50runs.txt' 

   # Random qrtr pixel offset for both
# fname = args.subdir + 'offsetBothRandomQrtrPixelLorR_deblendingTests_50runs.txt'

   # Random half pixel offset for both
# fname = args.subdir + 'offsetBothRandomHalfPixelLorR_deblendingTests_50runs.txt'
# fname = args.subdir + 'offsetBothRandomHalfPixelLorR_deblendingTests_10runs.txt'

   # Random qrtr pixel offset left, right one fixed
#fname = args.subdir + 'offsetLeftOneQrtrPixelLorRRightOneFixed_deblendingTests_50runs.txt'


   # Horiz sep- exact 0 sep and on points; and vert sep- exact 2" sep and on points
#fname = args.subdir + 'offsetEachVerticallyUpAndDownOneArcsecFromCenterAndNoHorizSep_deblendingTests_10runs.txt'

   # Horiz sep- exact 2" sep and on points; and vert sep- exact 2" sep and on points
#fname = args.subdir + 'offsetEachVerticallyUpAndDownOneArcsecFromCenterAndHorizSep2arcsec_deblendingTests_10runs.txt'

   # Vert sep
#fname =  args.subdir + 'offsetEachVerticallyUpAndDownRandomHalfPixelFromCenterAndHorizSep2arcsec_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetEachVerticallyUpRandomHalfPixelFromCenterAndHorizRandomHalfPixelFromCenter_deblendingTests_50runs.txt'

   # Horiz sep- exact 2" sep and on points; and vert displacement: +eps for both
#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0.001)__peak_B_(1.0, 0.001)_5_runs.txt'
#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0.001)__peak_B_(1.0, 0.001)_50_runs.txt'

  # Horiz sep- exact 2" sep and on points; and vert displacement: -eps for both
#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, -0.001)__peak_B_(1.0, -0.001)_5_runs.txt'
fname = args.subdir + 'deblendingTests_peak_A_(-1.0, -0.001)__peak_B_(1.0, -0.001)_50_runs.txt'

# fname = args.subdir + 'offsetVertSep2arcsecAndEachVerticallyUpAndDownRandomHalfPixelFromCenterAndHorizSep0arcsec_deblendingTests_50runs.txt'

   # Horiz sep- Random half pixel offset for both; and vert sep- also random half pixel offset for both


# fname = args.subdir + 'deblendingTests_peak_A_(-1, 0)__peak_B_(1, 0)_50_runsAndRandomOffsetHalfPixelEach.txt'
# fname = args.subdir + 'deblendingTests_peak_A_(0, -1)__peak_B_(0, 1)_50_runsAndRandomOffsetHalfPixelEach.txt'

#################################### Load up data
fitdat = np.loadtxt(fname)

fnumvec = fitdat[:,0]
fnum =  fnumvec

try:

    print "*********************** Has e2 in output file"
    # Initze all a vecs
    e1a_in = fitdat[:,1]
    e2a_in = fitdat[:,2]
    e1a_unbl = fitdat[:,3] 
    e1a_debl = fitdat[:,4] 
    e2a_unbl = fitdat[:,5] 
    e2a_debl = fitdat[:,6] 
    e1a_unblresid =  e1a_unbl - e1a_in 
    e1a_deblresid =  e1a_debl - e1a_in 
    
    e1b_in = fitdat[:,7]
    e2b_in = fitdat[:,8]
    e1b_unbl = fitdat[:,9] 
    e1b_debl = fitdat[:,10]
    e2b_unbl = fitdat[:,11] 
    e2b_debl = fitdat[:,12]
    e1b_unblresid =  e1b_unbl - e1b_in 
    e1b_deblresid =  e1b_debl - e1b_in 

except:

    print "*********************** No e2 in output file, reading just the orig files"
    e1a_in = fitdat[:,1]
    e2a_in = fitdat[:,2]
    e1a_unbl = fitdat[:,3] 
    e1a_debl = fitdat[:,4] 
    e1a_unblresid =  e1a_unbl - e1a_in 
    e1a_deblresid =  e1a_debl - e1a_in 

    e1b_in = fitdat[:,5]
    e2b_in = fitdat[:,6]
    e1b_unbl = fitdat[:,7] 
    e1b_debl = fitdat[:,8]
    e1b_unblresid =  e1b_unbl - e1b_in 
    e1b_deblresid =  e1b_debl - e1b_in 



# print 'fnum = ', fnum 

# Taking out the print statements for now
'''
print 'e1a_in = ' ,e1a_in 
print 'e1a_unbl  = ' ,e1a_unbl 
print 'e1a_debl = ', e1a_debl
print 'e1a_unblresid = ',e1a_unblresid

print 'e1b_in = ' ,e1b_in 
print 'e1b_unbl  = ' ,e1b_unbl 
print 'e1b_debl = ', e1b_debl
print 'e1b_unblresid = ',e1b_unblresid
'''

numfiles = 50

#e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

e1a_range = [0.5,  0,  -0.5]
e1b_range = [0.5,  0,  -0.5]

xshift = [0.01, 0.01, 0.01]

e1shifted = np.array(e1a_range) + np.array(xshift)
nbins = 10 # Set num of bins for histo
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

    # Run over all e1a fits
    for e1ain in e1a_range:
        e1astr = str(e1ain)
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
#        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
 #       print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        vece1a_in.append( e1a_in[i].mean() )
        vece1a_unbl.append( e1a_unbl[i].mean() )
        vece1a_debl.append( e1a_debl[i].mean() )
        vece1a_unblresid.append (e1a_unblresid[i].mean()  )
        vece1a_deblresid.append (e1a_deblresid[i].mean() )
        vece1a_unblerr.append( e1a_unbl[i].std()  )
        vece1a_deblerr.append( e1a_debl[i].std() ) 

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

#### Unblended fit plots
# (Note we are horizontally offsetting these points by xshift, as defined above)
    plt.scatter( e1shifted , vece1a_unblresid, color = 'g' , s=50.0  ) # The central value point -- 's' here is the size of the point -- 
    gline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1a" , linewidth=4.0)
    yline = plt.errorbar( e1shifted, vece1a_unblresid, vece1a_unblerr/np.sqrt(numfiles),  ecolor='y',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

#### Deblended fit plots
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
#        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
 #       print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        vece1b_in.append( e1b_in[i].mean() )
        vece1b_unbl.append( e1b_unbl[i].mean() )
        vece1b_debl.append( e1b_debl[i].mean() )
        vece1b_unblresid.append (e1b_unblresid[i].mean()  )
        vece1b_deblresid.append (e1b_deblresid[i].mean() )
        vece1b_unblerr.append( e1b_unbl[i].std()  )
        vece1b_deblerr.append( e1b_debl[i].std() ) 

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
