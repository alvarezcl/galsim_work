# Start: 10-24-2014
# MSSG

# Takes the files on disk and plots them up 

import matplotlib.pyplot as plt
import sys
import numpy as np
from argparse import ArgumentParser
import matplotlib.gridspec as gridspec

nonRoundObjs = False  # Boolean if we are plotting results of fits on more than a single type of obj (e.g. more than round ones)

########### Start the input arg parsing

parser = ArgumentParser()
parser.add_argument("--subdir", default="deblendsOutput/", help="output text filename")
args = parser.parse_args()

########## Flags

### Flag to tell if if we want to show the simfit columns (if the input file has them)
plot_simfit_flag = 1
multrow = 0 # If we want more than one row of plots -- for the "Slam" talks, set this to False; we'll use the same flag to remove the larger error bar, and only leave the error on the mean


############ Read in file
# fname = args.subdir + 'deblendingTests_50runs.txt' # Orig

#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_10_runs.txt'

#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_5_runs.txt'

#fname = args.subdir + 'offsetEachQuarterPixelAwayFromCenter_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBothQuarterPixelLeft_deblendingTests_50runs.txt'
#fname = args.subdir + 'offsetBoth0.005PixelLeft_deblendingTests_50runs.txt'     # 
#fname = args.subdir + 'offsetBoth0.005PixelRight_deblendingTests_50runs.txt'


 # Horiz sep- move centers to -eps and -eps
#fname = args.subdir + 'offsetBothOne200thPixelLeft_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to -eps and +eps
#fname = args.subdir + 'offsetLeftOne200thPixelLeftRightOne200thPixelRight_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to +eps and -eps
#fname = args.subdir + 'offsetLeftOne200thPixelRightRightOne200thPixelLeft_deblendingTests_50runs.txt' 

   # Horiz sep- move centers to +eps and +eps
#fname = args.subdir + 'offsetBothOne200thPixelRight_deblendingTests_50runs.txt' 

   # Random qrtr pixel offset for both
# fname = args.subdir + 'offsetBothRandomQrtrPixelLorR_deblendingTests_50runs.txt'

   # Random half pixel offset for both
# fname = args.subdir + 'offsetBothRandomHalfPixelLorR_deblendingTests_10runs.txt'


###################### Rounding vs. Truncation -- 3-24-2015
#fname = args.subdir + 'offsetBothRandomHalfPixelLorR_deblendingTests_50runs.txt'
#fname = args.subdir + 'roundedVsTruncated_offsetBothRandomHalfPixelLorR_deblendingTests_50runs.txt'
#fname = args.subdir + 'roundedVsTruncated_BothHalfpixelRandomShiftRight_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.txt'
#fname = args.subdir + 'roundedVsTruncated_BothQrtrpixelRandomShiftRight_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.txt'

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
#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, -0.001)__peak_B_(1.0, -0.001)_50_runs.txt'

# fname = args.subdir + 'offsetVertSep2arcsecAndEachVerticallyUpAndDownRandomHalfPixelFromCenterAndHorizSep0arcsec_deblendingTests_50runs.txt'

   # Horiz sep- Random half pixel offset for both; and vert sep- also random half pixel offset for both


#fname = args.subdir + 'deblendingTests_peak_A_(-1, 0)__peak_B_(1, 0)_50_runsAndRandomOffsetHalfPixelEach.txt'
#fname = args.subdir + 'deblendingTests_peak_A_(0, -1)__peak_B_(0, 1)_50_runsAndRandomOffsetHalfPixelEach.txt'


############################## Using arb point rotation (JEM new deblender, Jan 2015)
#fname = args.subdir + 'deblendingTests_peak_A_(-1.001, 0)__peak_B_(0.999, 0)_50_runs.InterpFlagTrue.txt'

fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.InterpFlagTrueANDForceInterpTrue.txt'

#fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.InterpFlagFalse.txt' # StandardCenters

########################### Round objs only (2/1/2015)
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_500_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_500_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_50_runs_withObjCenters.txt'  # 2-9-2015
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_500_runs_withObjCenters.txt'  # 2-9-2015
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_andRandomHalfPixelLorRoffset_5_runs_withObjCenters.txt'  # 2-10-2015
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_andRandomHalfPixelLorRoffset_500_runs_withObjCenters.txt'  # 2-10-2015

#fname = args.subdir + 'fdf_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50_runs.txt' # 4-3-2015

#fname = args.subdir + 'fdf_both.ellips.0.5.only.deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50_runs.txt'

#fname = args.subdir + 'fdf_run2_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50_runs.txt' # 4-3-2015

#fname =  'tmpdir/fdf_exactcenters_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50runs.txt' # 4-16-2015
#fname =  'tmpdir/fdf_fitcenters_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50runs.txt' # 4-16-2015

#fname =  'tmpdir/fdf_fitcenters_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_50runs.4-30-2015.txt' # 4-30-2015

fname =  'tmpdir/fdf_fitcenters_deblendingTests_peak_A_[-1  0]__peak_B_[1 0]_752runs.4-30-2015.txt' # 4-30-2015

# fname = args.subdir + 'deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.txt' # 7/7/2015

#fname =  'tmpdir/singletmpfile.txt'


#################################### Load up data
fitdat = np.loadtxt(fname)

fnumvec = fitdat[:,0]
fnum =  fnumvec

try:

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

    x0a_unbl = fitdat[:,13];    y0a_unbl = fitdat[:,14]
    x0a_debl = fitdat[:,15];    y0a_debl = fitdat[:,16]

    x0b_unbl = fitdat[:,17];    y0b_unbl = fitdat[:,18]
    x0b_debl = fitdat[:,19];    y0b_debl = fitdat[:,20]

    x0_a_sf = fitdat[:,21]
    x0_b_sf  = fitdat[:,22]
    e1a_sf  = fitdat[:,23]
    e1b_sf  = fitdat[:,24]
    e1a_sfresid =  e1a_sf - e1a_in 
    e1b_sfresid =  e1b_sf - e1b_in 

    peak_x0_a = fitdat[:,25]
    peak_y0_a  = fitdat[:,26]
    peak_x0_b = fitdat[:,27]
    peak_y0_b = fitdat[:,28]
    
    print "*********************** Has e2 in output file, and object centers"

except:   ############### NB: *check* the input file columns if there is a problem with the plots, they changed format from 8 to 12 columns at one point, so right now i put in the 12 column default

    print "*********************** No e2 in output file, reading just the orig files"
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

print 'fnum = ', fnum 


# Taking out the print statements for now
print 'e1a_in = ' ,e1a_in 
print 'e1a_unbl  = ' ,e1a_unbl 
print 'e1a_debl = ', e1a_debl
print 'e1a_unblresid = ',e1a_unblresid

print 'e1b_in = ' ,e1b_in 
print 'e1b_unbl  = ' ,e1b_unbl 
print 'e1b_debl = ', e1b_debl
print 'e1b_unblresid = ',e1b_unblresid


#sys.exit()

####################################### Initzn vars
#### How many runs we made
numfiles = 752
#numfiles = 50

#### Normal range i've been using
e1a_range = [0.5,  0, -0.5]
e1b_range = [0.5,  0, -0.5]

#### Extended range
#e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

#### To do just round ones, 2/1/2015
#e1a_range = [ 0]
#e1b_range = [ 0]

#### To separate data points and show error bars clearly
xshift = [0.01, 0.01, 0.01]
shiftfac = 3
e1shiftedL = np.array(e1a_range) - shiftfac *np.array(xshift)
e1shiftedR = np.array(e1a_range) + shiftfac *np.array(xshift)

#### Set num of bins for histo
nbins = 20 

###################### Initialize the geometry of the grid for the fig
gsvert = 20
gshoriz = 30
gs = gridspec.GridSpec(gsvert, gshoriz)
# For later when we put the plots in dif gridspec spots
gsh = gshoriz/3 
gsv = gsvert / 2
horizbuffer = 3
vertbuffer = 1

if nonRoundObjs: 
    totfig = plt.figure(figsize=(15,12))
figindex = 0

######### Print to screen (or ipy nb file)
print " Using file ", fname

print "\n\n **** About to plot fits to obj A now"

############################################################ e1 a plots
for e1bin in e1b_range:
    e1bstr = str(e1bin)

    vece1a_in = []  
    vece1a_unbl = []
    vece1a_debl = []
    vece1a_sf = []
    vece1a_unblresid = []  
    vece1a_deblresid = []  
    vece1a_sfresid = []  
    vece1a_unblerr = []  
    vece1a_deblerr = []  
    vece1a_sferr = []  


    # Run over all e1a fits
    for e1ain in e1a_range:
        e1astr = str(e1ain)
        # Pick the entries where the input vals are as below
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
        #print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
        #print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        vece1a_in.append( e1a_in[i].mean() )
        vece1a_unbl.append( e1a_unbl[i].mean() )
        vece1a_debl.append( e1a_debl[i].mean() )
        vece1a_unblresid.append (e1a_unblresid[i].mean()  )
        vece1a_deblresid.append (e1a_deblresid[i].mean() )
        vece1a_unblerr.append( e1a_unbl[i].std()  )
        vece1a_deblerr.append( e1a_debl[i].std() ) 

        if (plot_simfit_flag > 0):
            vece1a_sf.append( e1a_sf[i].mean() )
            vece1a_sfresid.append (e1a_sfresid[i].mean() )
            vece1a_sferr.append( e1a_sf[i].std() )         
            
    if nonRoundObjs:
        ########################################### e1 a plots

        ############ Plot it in this fig slot for this pass
        thisfig = totfig.add_subplot(gs[0:gsv-vertbuffer , figindex*gsh:(figindex+1)*gsh-horizbuffer])
#        xlimit = 0.55;    ylimit = 0.05
#        plt.xlim( -xlimit, xlimit);           plt.ylim( -ylimit, ylimit)

        print "  \n\n\n vece1a_unblresid =  ", vece1a_unblresid

        #### Unblended fit plots
        # (Note we are horizontally offsetting these points by xshift, as defined above)
        plt.scatter( e1a_range , vece1a_unblresid, color = 'y' , s=50.0  ) # The central value point -- 's' here is the size of the point 
        if multrow:
            gline = plt.errorbar( e1a_range, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid e1a for unblended objs" , linewidth=4.0)
        yline = plt.errorbar( e1a_range, vece1a_unblresid, vece1a_unblerr/np.sqrt(numfiles),  ecolor='y',linestyle=' ', label = r"Resid e1a for unblended objs, $\frac{\sigma}{ \sqrt{N} }$" , linewidth= 8.0 )

        #### Deblended fit plots
        plt.scatter( e1shiftedL, vece1a_deblresid, color = 'r' , s=50.0  )
        if multrow:
            bline = plt.errorbar(e1shiftedL , vece1a_deblresid, vece1a_deblerr,  ecolor='b',linestyle=' ', label = r"Resid e1a for deblended objs" , linewidth= 4.0)
        rline = plt.errorbar(e1shiftedL , vece1a_deblresid, vece1a_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = r"Resid e1a for deblended objs,  $\frac{\sigma}{ \sqrt{N} }$" , linewidth= 8.0 )

        #### Simfit plots
        if (plot_simfit_flag > 0):
            plt.scatter( e1shiftedR, vece1a_sfresid, color = 'c' , s=50.0  )
            if multrow:
                kline = plt.errorbar(e1shiftedR , vece1a_sfresid, vece1a_sferr,  ecolor='k',linestyle=' ', label = "Resid e1 for simfits" , linewidth= 4.0)
            cline = plt.errorbar(e1shiftedR , vece1a_sfresid, vece1a_sferr/np.sqrt(numfiles),  ecolor='c',linestyle=' ', label = r"Resid for simfit,   $\frac{\sigma}{ \sqrt{N} }$" , linewidth= 8.0 )

        # Title + legend
        plt.title("Resids for fits with e1bin = " +e1bstr )
        # plt.legend(loc=4,prop={'size':9}) # This loc is the lower right corner, and this is good font size for the box
        plt.xlabel('$e_{1a in}$',fontsize=18)
        plt.ylabel('$e_{1a fit}-e_{1a in}$',fontsize=18)
        plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
        
        # Increment figindex for next pass through this loop
        figindex += 1

#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")
plt.legend(loc=0,prop={'size':9}) # This loc=1 is the upper right corner, and this is good font size for the box
print "Finished e1a plots"

#################### Reset FigIndex for next loop
figindex = 0



############################################################ e1 b plots
for e1ain in e1a_range:
    e1astr = str(e1ain)

    vece1b_in = []  
    vece1b_unbl = []
    vece1b_debl = []
    vece1b_sf = []
    vece1b_unblresid = []  
    vece1b_deblresid = []  
    vece1b_sfresid = []  
    vece1b_unblerr = []  
    vece1b_deblerr = []  
    vece1b_sferr = []  


    for e1bin in e1b_range:
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
        #print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
        #print 'e1b_unbl, e1a_unbl = ', e1b_unbl[i], e1a_unbl[i] 

        vece1b_in.append( e1b_in[i].mean() )
        vece1b_unbl.append( e1b_unbl[i].mean() )
        vece1b_debl.append( e1b_debl[i].mean() )
        vece1b_unblresid.append (e1b_unblresid[i].mean()  )
        vece1b_deblresid.append (e1b_deblresid[i].mean() )
        vece1b_unblerr.append( e1b_unbl[i].std()  )
        vece1b_deblerr.append( e1b_debl[i].std() ) 

        if (plot_simfit_flag > 0):
            vece1b_sf.append( e1b_sf[i].mean() )
            vece1b_sfresid.append (e1b_sfresid[i].mean() )
            vece1b_sferr.append( e1b_sf[i].std() ) 
        
    if nonRoundObjs:
        ######################################################## e1 b plots
        print 'vece1b_unbl = ', vece1b_unbl
        #   plt.figure(figsize=(15,12))
        ############ Plot it in this fig slot for this pass

###### Show these only in the case of only one row of plots
        if (multrow):
            thisfig = totfig.add_subplot(gs[gsv:2*gsv-1 , figindex*gsh:(figindex+1)*gsh-horizbuffer])


            #### Unblended fit plots
            # (Note we are horizontally offsetting these points by xshift, as defined above)
# 7/19/2015
#            xlimit = 0.55;    ylimit = 0.05
#            plt.xlim( -xlimit, xlimit);            plt.ylim( -ylimit, ylimit)

            print "  \n\n\n vece1b_unblresid =  ", vece1b_unblresid

            plt.scatter( e1b_range , vece1b_unblresid, color = 'g' , s=50.0  )
            gline = plt.errorbar( e1b_range, vece1b_unblresid, vece1b_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1b" , linewidth=4.0)
            yline = plt.errorbar( e1b_range, vece1b_unblresid, vece1b_unblerr/np.sqrt(numfiles),  ecolor='y',linestyle=' ', label = "Resid for unblended fit, error/sqrt(N)" , linewidth= 8.0 )

            #### Deblended fit plots
            plt.scatter( e1shiftedL, vece1b_deblresid, color = 'b' , s=50.0  )
            bline = plt.errorbar(e1shiftedL, vece1b_deblresid, vece1b_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1b" , linewidth= 4.0)
            rline = plt.errorbar(e1shiftedL, vece1b_deblresid, vece1b_deblerr/np.sqrt(numfiles),  ecolor='r',linestyle=' ', label = "Resid for deblended fit, error/sqrt(N)" , linewidth= 8.0 )

            #### Simfit plots
            if (plot_simfit_flag > 0):
                plt.scatter( e1shiftedR, vece1b_sfresid, color = 'm' , s=50.0  )
                bline = plt.errorbar(e1shiftedR , vece1b_sfresid, vece1b_sferr,  ecolor='k',linestyle=' ', label = "Resid for simfit for e1b" , linewidth= 4.0)
                rline = plt.errorbar(e1shiftedR, vece1b_sfresid, vece1b_sferr/np.sqrt(numfiles),  ecolor='c',linestyle=' ', label = "Resid for simfit, error/sqrt(N)" , linewidth= 8.0 )

            plt.title("Resids for fits with e1ain = " +e1astr )
            plt.xlabel('$e_{1b in}$',fontsize=18)
            plt.ylabel('$e_{1b fit}-e_{1b in}$',fontsize=18)
            plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     


        '''
        #### Legend
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box
        '''

        figindex += 1
print "Finished e1b plots"

#totfig.legend()

plt.show()

sys.exit()
