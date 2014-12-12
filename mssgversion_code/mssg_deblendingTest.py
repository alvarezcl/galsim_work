# MSSG, based on JEM and LCA code
# 10/15/2014

import galsim
import numpy as np
import mssg_deblend
import mssg_drawLibrary
import lmfit
import ipdb ; import sys
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import random

################### Initze

#### Level to do printing at (setting it lower will print more stuff)
presetval = 1

bd=galsim.BaseDeviate(0) # Random num seed

############## Function to draw and return simple single gal img
def drawgal(peak =(0,0), e1 = 0, e2 = 0 , fwhm=1.0, flux=1.0e5,  psfshr=0, psfbeta=3.0 , psffwhm= 0.85):
    # Create gaussian gal objs, sheared in various directions
    gal = galsim.Gaussian(fwhm=fwhm, flux=flux).shear(g1=e1, g2= e2).shift(peak)
    # Add psf 
    psf = galsim.Moffat(beta=psfbeta, fwhm=psffwhm).shear(e1=psfshr,  e2= 0.0)
    convgal = galsim.Convolve([gal,psf])
    # Now make img
    pixelsize = 49;    pixelscale = 0.2
    # Img
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    img = convgal.drawImage(image=proto_image, method='phot', rng=bd)
    img.array[np.where(img.array < 0)] = 0.
    return img

############## Function to create the blended img

def create_blend(peak_a, peak_b, e1a = 0, e1b = 0 , e2a = 0, e2b = 0):
    # Create gaussian gal objs, sheared in various directions

    gal1 = galsim.Gaussian(fwhm=1.0, flux=100000.0).shear(g1=e1a, g2= e2a).shift(peak_a)
    gal2 = galsim.Gaussian(fwhm=1.0, flux=100000.0).shear(g1=e1b, g2= e2b).shift(peak_b)
    
    # Add psf 
    psfshr = 0.00
    psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1=psfshr,  e2=-0.0)

    convgal1 = galsim.Convolve([gal1,psf])
    convgal2 = galsim.Convolve([gal2,psf])
    
#    convgal1 = gal1 
#    convgal2 = gal2

    # Now make imgs
    pixelsize = 49
    pixelscale = 0.2

    # Img1
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    #image1 = gal1.drawImage(image=proto_image, method='fft')
    #image1 = gal1.drawImage(image=proto_image, method='phot', rng=bd)
    # With psf convln
    #image1 = convgal1.drawImage(image=proto_image, method='fft')
    image1 = convgal1.drawImage(image=proto_image, method='phot', rng=bd)
    image1.array[np.where(image1.array < 0)] = 0.



    if plotflag > presetval:
    #        print " >>>>>> Plotting img1, first pass"
#        img1 =        drawgal(peak = (-1,0), e1 = e1a, e2 = e2a, psfbeta = 3.0 , psffwhm = 0.85, flux = 1.0e5, fwhm=1.0 )
 #       plt.imshow( img1.array , origin='lower');    
  #      plt.colorbar()
   #     plt.show()
        plt.title(" FIRST PLOT:  Img obj a")
        print " >>>>>> Plotting img obj a"
        plt.imshow( image1.array , origin='lower');    
        plt.colorbar()
        plt.show()

#        print " img1.array.sum() " , img1.array.sum()
        print " image1.array.sum() " , image1.array.sum()

###        print " >>>>>> Plotting img1 - image1"
   #     plt.imshow( img1.array - image1.array , origin='lower');    
    #    plt.colorbar()
     #   plt.show()

    # Img2
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    #    image2 = gal2.drawImage(image=proto_image, method='fft')
    #image2 = gal2.drawImage(image=proto_image, method='phot', rng=bd)
    # With psf convln
    #image2 = convgal2.drawImage(image=proto_image, method='fft')
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=bd)
    image2.array[np.where(image2.array < 0)] = 0.


    
    if plotflag > presetval: 
        plt.title(" Img obj b")
        print " >>>>>> Plotting img obj b"
        plt.imshow( image2.array , origin='lower');    
        plt.colorbar()
        plt.show()

### Add them into one image
    imagesum =  image1+image2

    if plotflag > presetval:
        plt.title(" Img obj a + b")
        print " >>>>>> Plotting imgsum"
        plt.imshow( imagesum.array , origin='lower');    
        plt.colorbar()
        plt.show()

    return imagesum, [image1, image2]







############################################################################# Main
if __name__ == '__main__':


# Parse command line args
    parser = ArgumentParser()
    parser.add_argument("--outfile", default="deblendsOutput/deblendingTests", help="output text filename")
    parser.add_argument("--e1a", default=0, type=float, help="e1a in")
    parser.add_argument("--e2a", default=0, type=float, help="e2a in")
    parser.add_argument("--e1b", default=0, type=float, help="e1b in")
    parser.add_argument("--e2b", default=0, type=float, help="e2b in")
    parser.add_argument("--plotflag", default=0, type=int, help="Set to 1 to make plots")

    args = parser.parse_args()

# Get ellips from cmd line
    e1a_in = args.e1a
    e2a_in = args.e2a

    e1b_in = args.e1b
    e2b_in = args.e2b

    plotflag = args.plotflag

# Centroids
#    peak_a = (-1.001,0);   peak_b = (0.999,0)  # Horiz sep- move centers to -eps and -eps
#    peak_a = (-1.001,0);   peak_b = (1.001,0)  # Horiz sep- move centers to -eps and +eps
#    peak_a = (-0.999,0);   peak_b = (1.001,0)  # Horiz sep- move centers to +eps and +eps
#    peak_a = (-0.999,0);   peak_b = (0.999,0)  # Horiz sep- move centers to +eps and -eps

#    peak_a = (-1.05,0);   peak_b = (0.95,0)  # Horiz sep- move centers both to left quarter pixel
#    peak_a = (-1.05,0);   peak_b = (1.05,0)  # Horiz sep - move centers apart quarter pixel
#    peak_a = (-1.0,0);   peak_b = (1.0,0)    # Horiz sep - centers separated by 2", EXACTLY

#    peak_a = (0,1);   peak_b = (0,-1)  # Vert sep

# Convert peak to pixels
#    peak_a = (-0.999,0);   peak_b = (0.999,0)  # Horiz sep- move centers to +eps and -eps

# Vertical shifts
    peak_a = (0,-1.0);   peak_b = (0,1.0)    # Horiz sep - 0, and vertical shift by 2.0"
#   peak_a = (-1.0,-1.0);   peak_b = (1.0,1.0)    # Horiz sep - centers separated by 2", EXACTLY, and vertical shift by 2.0"


####### Convert to pixels
#    peaks_pix = [[p1/0.2 for p1 in peak_a],
#                 [p2/0.2 for p2 in peak_b]]

#    print peaks_pix 

################################################## Start loops

    fitdat = []

#    e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#    e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

    e1a_range = [0.5,  0, -0.5]
    e1b_range = [0.5,  0, -0.5]

    e2ain = 0
    e2bin = 0

    xavec =[] ;     xbvec =[]


# Set size of partial pixels
    qrtrpixel = 0.25
    halfpixel = 0.5

    print " \n\n\n peak_a = ",  peak_a 

    '''
# Random number test loop   
    for filenum in xrange(0,numfiles):
        
        xashift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist
        xbshift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist
#        print "xashift = ", xashift
        
        xashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
        xbshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist

        xavec.append(xashift) ;         xbvec.append(xbshift)
        
        xav = np.array(xavec) ;     xbv = np.array(xbvec)

    print " mean(xashift) = ", xav.mean() 
    print " std(xashift) = ", xav.std() 
    '''
     
#    sys.exit()
    
    numfiles = 50  # Number of runs
    for filenum in xrange(0,numfiles):

### Run over ellips
        for e1bin in e1b_range:
            for e1ain in e1a_range:

                print " ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin

#                xashift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist
#                xbshift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist

                xashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
                xbshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist

                yashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
                ybshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist

#                xbshift = 0
                xaoff = xashift/5 ; xboff = xbshift/5  # Converts to arcsec, by div by 5
                yaoff = yashift/5 ; yboff = ybshift/5  # Converts to arcsec, by div by 5
#                peak_a = (-1.0 + xaoff,-1.0 + yaoff);   peak_b = (1.0 + xboff,1.0 + yboff)  # Horiz sep- move centers by the random offsets, in arcsec,  Vert sep- also move centers by the random offsets, in arcsec
#                peak_a = (-1.0 , yaoff);   peak_b = (1.0 , yboff)  # Horiz sep- move centers by the random offsets, in arcsec,  Vert sep- also move centers by the random offsets, in arcsec
                peak_a = (0 , -1.0 + yaoff);   peak_b = (0 , 1.0 + yboff)  # Horiz sep- move centers by the random offsets, in arcsec,  Vert sep- also move centers by the random offsets, in arcsec
                peaks_pix = [[p1/0.2 for p1 in peak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in peak_b]]
                
                print " Arcsec: peaks = " , peak_a, peak_b 
                print " Pixels: peaks_pix = " ,  peaks_pix 

# Create the blended object using funct above 
#   blend = single image where we've added both
#   unblends = vector of 2 imgs, each of the ind ones
                blend, unblends = create_blend(peak_a, peak_b, e1a = e1ain,  e2a = e2ain, e1b = e1bin ,e2b = e2bin)

                if plotflag > presetval:
                    plt.title(" Img blended obj - (a+b) ")
                    print " >>>>>> Plotting blend.array - (unblends[0].array + unblends[1].array)  "
                    plt.imshow( blend.array - (unblends[0].array + unblends[1].array) , origin='lower');    
                    plt.colorbar()
                    plt.show()

# Use deblending code to separate them
#    templates = for each img
#    template_fractions
#    children = vector of 2 imgs, best estimates from deblending code
                templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix)

                if plotflag > presetval:
                    plt.title(" Template obj a ")
                    print " >>>>>> Plotting template fraction a"
                    plt.title(" Template 1")
                    plt.imshow( template_fractions[0] , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Template obj b ")
                    print " >>>>>> Plotting template fraction b"
                    plt.imshow( template_fractions[1] , origin='lower');    
                    plt.colorbar()
                    plt.show()

    ######### Plot children
                if plotflag > presetval:
                    print " >>>>>> Plotting Unblended img a "
                    plt.title(" Unblended img a ")
                    plt.imshow( unblends[0].array , origin='lower');            plt.colorbar();        plt.show()
                    print " >>>>>> Plotting Deblended child a "
                    plt.title(" Deblended child a ")
                    plt.imshow( children[0] , origin='lower');            plt.colorbar();        plt.show()
                    print " >>>>>> Plotting Resid of: (Deblended child a - Unblended img a) "
                    plt.title("Resid of: (Deblended child a - Unblended img a)  ")
                    plt.imshow( children[0] - unblends[0].array , origin='lower');            plt.colorbar();        plt.show()

                    print " >>>>>> Plotting Unblended img b "
                    plt.title(" Unblended img b ")
                    plt.imshow( unblends[1].array , origin='lower');            plt.colorbar();        plt.show()
                    print " >>>>>> Plotting Deblended child b "
                    plt.title(" Deblended child b ")
                    plt.imshow( children[1] , origin='lower');            plt.colorbar();        plt.show()
                    print " >>>>>> Plotting Resid of: (Deblended child b - Unblended img b) "
                    plt.title("Resid of: (Deblended child b - Unblended img b)  ")
                    plt.imshow( children[1] - unblends[1].array , origin='lower');            plt.colorbar();        plt.show()



    # Add the flux fraction imgs by hand 
    #   templatesum = adding two template_fraction imgs
                templatefracsum = template_fractions[0] + template_fractions[1]
    #    tarray = templatesum.array

                if plotflag > presetval:
                    plt.title(" Template fraction sum of a+b ")
                    print " >>>>>> Plotting templatesum"
                    plt.imshow( templatefracsum , origin='lower');    
                    plt.colorbar()
                    plt.show()


    # Add the child imgs by hand 
                sum_children = children[0] + children[1]    

                if plotflag > presetval:
                    plt.title(" Deblended children a + b")
                    print " >>>>>> Plotting sum_children "
                    plt.imshow( sum_children  , origin='lower');    
                    plt.colorbar()
                    plt.show()

        # Get resid
                residual = blend.array - sum_children

                if plotflag > presetval:
                    plt.title(" Residual: blended img - sum of deblended children ")
                    print " >>>>>> Plotting resid = blend.array - sum_children "
                    plt.imshow( residual  , origin='lower');    
                    plt.colorbar()
                    plt.show()

                print 'blend.array.sum() = ' , blend.array.sum()
                print 'children[0].sum() + children[1].sum() = ', children[0].sum() + children[1].sum()


    ######################## Run fits
                print " ***** Now about to run lmfit "
        # Common params to all
                imsize = 49
                pixel_scale = 0.2
                galtype = galsim.Gaussian

        # The below are just initial guesses for the fitter
                minmaxval = 0.7 # Temporary limit for the fitter till i fix it better by doing value in quadrature
                fit_params = lmfit.Parameters()  
                fit_params.add('e1_a', value=0.0, min= -minmaxval, max=minmaxval)
                fit_params.add('e2_a', value=0.0, min= -minmaxval, max=minmaxval)
                fit_params.add('flux_a', value=2000.0)
                fit_params.add('x0_a', value=0)
                fit_params.add('y0_a', value=0)
                fit_params.add('hlr_a', value=0.43)

    ####### Whether to convolve with PSF
                dopsfconvln = 'y'

    ##### Obj a
         ## Unbl Obj a        
                origimg = unblends[0]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype , dopsfconvln) ) 

        # Extract vals
                e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                print "\n\n ********************* Unbl Obj a "
                print "e1_a = ", e1_a, " ,  e1err = ", e1err
                print "e2_a = ", e2_a, " ,  e2err = ", e2err

                e1a_unbl = e1_a ; e2a_unbl = e2_a 
                e1a_unblerr = e1err ; e2a_unblerr = e2err 

    ################# Unbl obj a plots
                if plotflag > (presetval-1):
                    plt.title(" Unblended fit results for a")
                    print " >>>>>> Plotting result of unbl fit obj a "
                    fitobj_a = drawgal(peak = peak_a , e1 = e1_a, e2 = e2_a)
                    plt.imshow( fitobj_a.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj a")
                    print " >>>>>> Plotting orig obj a "
                    plt.imshow( origimg.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (unblended fit result for a - orig obj a )")
                    print " >>>>>> Plotting resid of (unbl fit obj a - orig obj a ) "
                    plt.imshow(  fitobj_a.array - origimg.array   , origin='lower');    
                    plt.colorbar()
                    plt.show()


        # Report the parameters to the interpreter screen                        
        #    lmfit.report_errors(mlresult.params)

    ############ Deblended Obj a
                origimg = children[0]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) ) 

        # Extract vals
                e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                print "\n *********  Deblended Obj a "
                print "e1_a = ", e1_a, " ,  e1err = ", e1err
                print "e2_a = ", e2_a, " ,  e2err = ", e2err

                e1a_debl = e1_a ; e2a_debl = e2_a 
                e1a_deblerr = e1err ; e2a_deblerr = e2err 

    # Report the parameters to the interpreter screen                        
    #    lmfit.report_errors(mlresult.params)

    ################# Debl obj a plots
                if plotflag > (presetval-1):
                    plt.title(" Deblended fit results for a")
                    print " >>>>>> Plotting result of debl fit obj a "
                    fitobj_a = drawgal(peak = peak_a , e1 = e1_a, e2 = e2_a)
                    plt.imshow( fitobj_a.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj a")
                    print " >>>>>> Plotting orig obj a "
                    plt.imshow( origimg  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (deblended fit result for a - orig obj a )")
                    print " >>>>>> Plotting resid of (orig obj a - debl fit obj a) "
                    plt.imshow( origimg - fitobj_a.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()

    ########### Unbl Obj b        
                origimg = unblends[1]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) )  

        # Extract vals
                e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                print "\n\n *********************  Unbl Obj b "
                print "e1_b = ", e1_b, " ,  e1err = ", e1err
                print "e2_b = ", e2_b, " ,  e2err = ", e2err

                e1b_unbl = e1_b ; e2b_unbl = e2_b 
                e1b_unblerr = e1err ; e2b_unblerr = e2err 

    ################# Unbl obj b plots
                if plotflag > (presetval-1):
                    plt.title(" Unblended fit results for b")
                    print " >>>>>> Plotting result of unbl fit obj b "
                    fitobj_b = drawgal(peak = peak_b , e1 = e1_b, e2 = e2_b)
                    plt.imshow( fitobj_b.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj b")
                    print " >>>>>> Plotting orig obj b "
                    plt.imshow( origimg.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (unblended fit result for b - orig obj b )")
                    print " >>>>>> Plotting resid of (orig obj b - unbl fit obj b) "
                    plt.imshow( origimg.array - fitobj_b.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()


        # Report the parameters to the interpreter screen                        
        #    lmfit.report_errors(mlresult.params)

    ######## Deblended Obj b
                origimg = children[1]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) )  

        # Extract vals  -- note that calling below: mlresult.params['e1_a']  is correct because the '_a' name is defined that way inside the fitter
                e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                print "\n **********  Deblended Obj b "
                print "e1_b = ", e1_b, " ,  e1err = ", e1err
                print "e2_b = ", e2_b, " ,  e2err = ", e2err

                e1b_debl  = e1_b ; e2b_debl  = e2_b 
                e1b_deblerr = e1err ; e2b_deblerr = e2err 

    ################# Debl obj b plots
                if plotflag > (presetval-1):
                    plt.title(" Deblended fit results for b")
                    print " >>>>>> Plotting result of unbl fit obj b "
                    fitobj_b = drawgal(peak = peak_b , e1 = e1_b, e2 = e2_b)
                    plt.imshow( fitobj_b.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj b")
                    print " >>>>>> Plotting orig obj b "
                    plt.imshow( origimg  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (deblended fit result for b - orig obj b )")
                    print " >>>>>> Plotting resid of (orig obj b - unbl fit obj b) "
                    plt.imshow( origimg - fitobj_b.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()

    ################### Result vec for this fit
                fitresults = [int(filenum), e1ain, e2ain, e1a_unbl, e1a_debl,   e1bin, e2bin, e1b_unbl, e1b_debl]
                fitdat.append(fitresults)
                print 'len(fitdat) = ', len(fitdat)

    ################################## End of all loops
    fitarray = np.array(fitdat)

    print(fitarray)

    np.savetxt('deblendingTestsvert_50runs.txt', fitarray, header="filenum   e1a_in e2a_in   e1a_unbl   e1a_debl    e1b_in e2b_in  e1b_unbl  e1b_debl")


