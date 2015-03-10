# Start: 10-24-2014
# MSSG

# Takes the files on disk and plots them up 

import matplotlib.pyplot as plt
import sys
import numpy as np
from argparse import ArgumentParser
import matplotlib.gridspec as gridspec
import triangle

nonRoundObjs = False  # Boolean if we are plotting results of fits on more than a single type of obj (e.g. more than round ones)

########### Start the input arg parsing

parser = ArgumentParser()
parser.add_argument("--subdir", default="deblendsOutput/", help="output text filename")
args = parser.parse_args()

############ Read in file
# fname = args.subdir + 'deblendingTests_50runs.txt' # Orig

########################### Round objs only (2/1/2015)
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_50_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-1.0, 0)__peak_B_(1.0, 0)_500_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_500_runs.txt'
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_50_runs_withObjCenters.txt'  # 2-9-2015
fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_500_runs_withObjCenters.txt'  # 2-9-2015
#fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_andRandomHalfPixelLorRoffset_5_runs_withObjCenters.txt'  # 2-10-2015
fname = args.subdir + 'roundObjsOnly_deblendingTests_peak_A_(-0.8, 0)__peak_B_(0.8, 0)_andRandomHalfPixelLorRoffset_500_runs_withObjCenters.txt'  # 2-10-2015


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

    print "*********************** Has e2 in output file, and object centers"

except:

    print "*********************** No e2 in output file -- exiting "
    sys.exit()



################################ To do just round ones, 2/1/2015
e1a_range = [ 0]
e1b_range = [ 0]

xshift = [0.01, 0.01, 0.01]

e1shifted = np.array(e1a_range) + np.array(xshift)
nbins = 20 # Set num of bins for histo

################################### Initialize the geometry of the grid for the fig
gs = gridspec.GridSpec(2,3)

################################### Print to screen (or ipy nb file)
print " Using file ", fname

print "\n\n **** About to plot fits to obj A now"

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

############### Declare the overall fig box
	totfig = plt.figure(figsize=(20,12))

        #  Make histo of e1a fit vals
        avg_e1a_unbl = round(np.mean(e1a_unbl),4)
        avg_e1a_debl = round(np.mean(e1a_debl),4) 
        thisfig = totfig.add_subplot(gs[0,0])
        plt.title("e1a unbl & debl fit for e1ain = " +e1astr + ", e1bin = " +e1bstr )
        xleft = -0.06
        xwidth = 0.12
        xright = xleft + xwidth
        plt.hist(e1a_unbl[i],range = [xleft,xright], bins= nbins,  alpha=0.5, label = "e1a unblended fit; avg = " + str(avg_e1a_unbl) )
        plt.hist(e1a_debl[i],range = [xleft,xright], bins= nbins,   alpha=0.5, label = "e1a deblended fit; avg = " + str(avg_e1a_debl) , color = 'g')
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box     
        print "Mean of e1a_unbl ",  avg_e1a_unbl
        print "Mean of e1a_debl ",  avg_e1a_debl

#        plt.show()

############## Histos of dists
        #  Make histo of x0a fit vals
        avg_x0a_unbl = round(np.mean(x0a_unbl),4)
        avg_x0a_debl = round(np.mean(x0a_debl),4) 
        thisfig = totfig.add_subplot(gs[0,1])
        plt.title("x0a unbl & debl fit for x0ain = -0.8 " )
        xleft = -0.95 
        xwidth = 0.25
        xright = xleft + xwidth
        plt.hist(x0a_unbl[i],range = [xleft,xright], bins= nbins,  alpha=0.5, label = "x0a unblended fit; avg = " + str(avg_x0a_unbl) )
        plt.hist(x0a_debl[i],range = [xleft,xright], bins= nbins, color = 'g',  alpha=0.5 , label = "x0a deblended fit; avg = " + str(avg_x0a_debl) )
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box
	print "Mean of x0a_unbl ", np.mean(x0a_unbl)
	print "Mean of x0a_debl ", np.mean(x0a_debl)
#        plt.show()

        #  Make histo of y0a fit vals
        avg_y0a_unbl = round(np.mean(y0a_unbl),4)
        avg_y0a_debl = round(np.mean(y0a_debl),4) 
        thisfig = totfig.add_subplot(gs[0,2])
        plt.title("y0a unbl  & debl fit dist for y0ain = 0 " )
        yleft = -0.008 
        ywidth = 0.016
        yright = yleft + ywidth
        plt.hist(y0a_unbl[i],range = [yleft,yright],bins= nbins,  alpha=0.5 , label = "y0a unblended fit; avg = " + str(avg_y0a_unbl) )
        plt.hist(y0a_debl[i],range = [yleft,yright],bins= nbins, color = 'g',  alpha=0.5 , label = "y0a deblended fit; avg = " + str(avg_y0a_debl) )
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box
	print "\n Mean of y0a_unbl ", np.mean(y0a_unbl)
	print "Mean of y0a_debl ", np.mean(y0a_debl)

print "\n\n **** About to plot fits to obj B now"

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


############## Histos of dists        
        #  Make histo of e1b fit vals
        avg_e1b_unbl = round(np.mean(e1b_unbl),4)
        avg_e1b_debl = round(np.mean(e1b_debl),4) 
        thisfig = totfig.add_subplot(gs[1,0])
        plt.title("e1b unbl & debl fit dist for e1ain = " +e1astr + ", e1bin = " +e1bstr )
        xleft = -0.06
        xwidth = 0.12
        xright = xleft + xwidth
        plt.hist(e1b_unbl[i],range = [xleft,xright], bins= nbins,  alpha=0.5, label = "e1b unblended fit; avg = " + str(avg_e1b_unbl) )
        plt.hist(e1b_debl[i], range = [xleft,xright], bins= nbins,   alpha=0.5, label = "e1b deblended fit; avg = " + str(avg_e1b_debl) , color = 'g')
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box     
        print "Mean of e1b_unbl ",  avg_e1b_unbl
        print "Mean of e1b_debl ",  avg_e1b_debl

#        plt.show()

        #  Make histo of x0b fit vals
        xright = +0.95
        xwidth = 0.25
        xleft = xright - xwidth

        avg_x0b_unbl = round(np.mean(x0b_unbl),4)
        avg_x0b_debl = round(np.mean(x0b_debl),4) 
        thisfig = totfig.add_subplot(gs[1,1])
        plt.title("x0a unbl  & debl fit dist for x0ain = -0.8 " )
        plt.hist(x0b_unbl[i],range = [xleft,xright], bins= nbins,  alpha=0.5, label = "x0a unblended fit; avg = " + str(avg_x0b_unbl) )
        plt.hist(x0b_debl[i],range = [xleft,xright], bins= nbins, color = 'g',  alpha=0.5 , label = "x0a deblended fit; avg = " + str(avg_x0b_debl) )
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box
	print "\n Mean of x0b_unbl ", np.mean(x0b_unbl)
	print "Mean of x0b_debl ", np.mean(x0b_debl)
#        plt.show()

        #  Make histo of y0b fit vals
        avg_y0b_unbl = round(np.mean(y0b_unbl),4)
        avg_y0b_debl = round(np.mean(y0b_debl),4) 
        thisfig = totfig.add_subplot(gs[1,2])
        plt.title("y0a unbl  & debl fit dist for y0ain = 0 " )
        plt.hist(y0b_unbl[i],bins= nbins,range = [yleft,yright],  alpha=0.5 , label = "y0a unblended fit; avg = " + str(avg_y0b_unbl) )
        plt.hist(y0b_debl[i],bins= nbins,range = [yleft,yright], color = 'g',  alpha=0.5 , label = "y0a deblended fit; avg = " + str(avg_y0b_debl) )
        plt.legend(loc=1,prop={'size':9}) # This loc is the upper right corner, and this is good font size for the box
        avg_y0b_unbl = np.mean(y0b_unbl)
        avg_y0b_debl = np.mean(y0b_debl)
	print "\n Mean of y0b_unbl ", np.mean(y0b_unbl)
	print "Mean of y0b_debl ", np.mean(y0b_debl)

# Show all the histos
        plt.show()

        
        y = [x0a_debl.tolist(),x0b_debl.tolist(), y0a_debl.tolist(), y0b_debl.tolist() ]
        z = np.array(y)
        fig = triangle.corner(z.T)
        plt.show()

