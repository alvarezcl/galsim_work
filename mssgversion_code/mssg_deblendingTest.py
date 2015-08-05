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
import triangle

################### Initze

#### Level to do printing at (setting it lower will print more stuff)
presetval = 1

bd=galsim.BaseDeviate(3) # Random num seed -- when set to zero, uses machine time

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

    ############# Img1
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    # With psf convln
    #image1 = convgal1.drawImage(image=proto_image, method='fft')  # If using fft method
    image1 = convgal1.drawImage(image=proto_image, method='phot', rng=bd)
    image1.array[np.where(image1.array < 0)] = 0.

    if plotflag > presetval:
        plt.title(" FIRST PLOT:  Img obj a")
        if v>0:
            print " >>>>>> Plotting img obj a"
        plt.imshow( image1.array , origin='lower');    
        plt.colorbar()
        plt.show()
        print " image1.array.sum() " , image1.array.sum()

    ############# Img2
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    # With psf convln
    #image2 = convgal2.drawImage(image=proto_image, method='fft')  # If using fft method
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=bd)
    image2.array[np.where(image2.array < 0)] = 0.
    
    if plotflag > presetval: 
        plt.title(" Img obj b")
        print " >>>>>> Plotting img obj b"
        plt.imshow( image2.array , origin='lower');    
        plt.colorbar()
        plt.show()

    ########### Add them into one image
    imagesum =  image1+image2

    if plotflag > presetval:
        plt.title(" Img obj a + b")
        print " >>>>>> Plotting imgsum"
        plt.imshow( imagesum.array , origin='lower');    
        plt.colorbar()
        plt.show()

    return imagesum, [image1, image2]




########################################################################################################### Main
if __name__ == '__main__':


# Parse command line args
    parser = ArgumentParser()
    parser.add_argument("--outdir", default="deblendsOutput/", help="Output text filename prepend")
    parser.add_argument("--e1a", default=0, type=float, help="e1a in")
    parser.add_argument("--e2a", default=0, type=float, help="e2a in")
    parser.add_argument("--e1b", default=0, type=float, help="e1b in")
    parser.add_argument("--e2b", default=0, type=float, help="e2b in")
    parser.add_argument("--verbosity", default=1, type=int, help="Set up to 3 to print max stuff to screen")
    parser.add_argument("--plotflag", default=0, type=int, help="Set to 1 to make plots")
    parser.add_argument("--numfile", default=1, type=int, help="Number of file configs to run and fit")
    parser.add_argument("--sepstep", default=0, type=float, help="Separation step in the pixel from L to R (should go 0->1)")

    args = parser.parse_args()

# Get ellips from cmd line
    e1a_in = args.e1a
    e2a_in = args.e2a

    e1b_in = args.e1b
    e2b_in = args.e2b

    plotflag = args.plotflag
    v = args.verbosity
    numfiles = args.numfile
    sepstep = args.sepstep
 
    outdir = args.outdir

######### Horiz shifts
    peak_a = (-1.0,0);   peak_b = (1.0 + sepstep/5 ,0)    # Horiz sep - centers separated by 2", EXACTLY, unless we've called this with a nonzero sepstep moving peak_b to the right (we only want to move across the pixel) 

#    peak_a = (-0.8,0);   peak_b = (0.8,0)    # Horiz sep - centers separated by 1.6", to match with Luis, EXACTLY
#    peak_a = (-1.001,0);   peak_b = (0.999,0)  # Horiz sep- move centers to -eps and -eps
#    peak_a = (-1.001,0);   peak_b = (1.001,0)  # Horiz sep- move centers to -eps and +eps
#    peak_a = (-0.999,0);   peak_b = (1.001,0)  # Horiz sep- move centers to +eps and +eps
#    peak_a = (-0.999,0);   peak_b = (0.999,0)  # Horiz sep- move centers to +eps and -eps

#    peak_a = (-1.05,0);   peak_b = (0.95,0)  # Horiz sep- move centers both to left quarter pixel
#    peak_a = (-1.05,0);   peak_b = (1.05,0)  # Horiz sep - move centers apart quarter pixel
#    peak_a = (-1.0,0);   peak_b = (1.0,0)    # Horiz sep - centers separated by 2", EXACTLY



######### Vertical shifts
#    peak_a = (0, -1.0);   peak_b = (0, 1.0)    # Horiz sep - 0, and vertical shift by 2.0"
#   peak_a = (-1.0,-1.0);   peak_b = (1.0,1.0)    # Horiz sep - centers separated by 2", EXACTLY, and vertical shift by 2.0"
#    peak_a = (-1.0,+0.001);   peak_b = (1.0,+0.001)    # Horiz sep - centers separated by 2", EXACTLY, and vertical shift up on each by +eps
#    peak_a = (-1.0,-0.001);   peak_b = (1.0,-0.001)    # Horiz sep - centers separated by 2", EXACTLY, and vertical shift up on each by -eps


####### Convert to pixels
#    peaks_pix = [[p1/0.2 for p1 in peak_a],
#                 [p2/0.2 for p2 in peak_b]]

#    print peaks_pix 

################################################## Start loops

    ### Initze the output var
    fitdat = []

#    e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
#    e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

#### Normal range i've been using
#    e1a_range = [0.5,  0, -0.5]
#    e1b_range = [0.5,  0, -0.5]

#### To do just round ones, 2/1/2015
    e1a_range = [ 0]
    e1b_range = [ 0]

    e2ain = 0
    e2bin = 0

    xavec =[] ;     xbvec =[]

# Set size of partial pixels
    qrtrpixel = 0.25
    halfpixel = 0.5

    if v>0:
        print " \n\n\n peak_a = ",  peak_a 

    '''
    # Random number test loop   
    numfiles = 50  # Number of runs
    for filenum in xrange(0,numfiles):
        
        xashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
        xbshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
        print "xashift = ", xashift

        xavec.append(xashift) ;         xbvec.append(xbshift)
        
        xav = np.array(xavec) ;     xbv = np.array(xbvec)

    print " mean(xashift) = ", xav.mean() 
    print " std(xashift) = ", xav.std() 
  
    sys.exit()
    '''
    

    origpeak_a = peak_a ; origpeak_b = peak_b
# origpeak_a = (-1, 0) ; origpeak_b = (1, 0)  # horiz offset, A is L, and B is R
#    origpeak_a = (0, -1) ; origpeak_b = (0, 1)  # vert offset, A is below, and B above


    for filenum in xrange(0,numfiles):

### Run over ellips
        for e1bin in e1b_range:
            for e1ain in e1a_range:
                
                if v>0:
                    print " ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin


### Initze all shifts to zero
                xashift = 0 ;                 xbshift = 0 
                yashift = 0 ;                 ybshift = 0 

#### For random offsets of the centroid, comment the following in
                ## Qrtr pixel offsets
                #                xashift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist
                #                xbshift = (random.random() / 2 ) - qrtrpixel # -0.25 to 0.25 flat dist

                ## Half pixel offsets - horiz
                # xashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
                # xbshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist

                ## Half pixel offsets just to the right
                # xashift = (random.random() / 2 )  # 0 to 0.5 flat dist
                # xbshift = (random.random() / 2 )  # 0 to 0.5 flat dist

                ## Qrtr pixel offsets just to the right
                # xashift = (random.random() / 4 )  # 0 to 0.25 flat dist
                # xbshift = (random.random() / 4 )  #  0 to 0.25 flat dist 

                ## Half pixel offsets - vert
#                yashift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist
 #               ybshift = (random.random() ) - halfpixel # -0.5 to 0.5 flat dist


                xaoff = xashift/5 ; xboff = xbshift/5  # Converts to arcsec, by div by 5
                yaoff = yashift/5 ; yboff = ybshift/5  # Converts to arcsec, by div by 5

                xashift = (xaoff, 0) ;                xbshift = (xboff, 0) # horizshift
                yashift = (0 , yaoff) ;                ybshift = (0 , yboff) # vertshift

#### End of  random offsets section

                if v>1:
                     print 'xashift = ', xashift
                     print 'xbshift = ', xbshift

                if v>1:
                    print  " np.array(origpeak_a) ,  np.array(xashift)  np.array(origpeak_a) + np.array(xashift) = ", np.array(origpeak_a) ,  np.array(xashift) , np.array(origpeak_a) + np.array(xashift)

# Horiz sep- move centers by the random offsets, in arcsec
                peak_a =  np.array(origpeak_a) + np.array(xashift);                peak_b =  np.array(origpeak_b) + np.array(xbshift )  

#  Vert sep- move centers by the random offsets, in arcsec
#                peak_a =  np.array(origpeak_a) + np.array(yashift);                peak_b =  np.array(origpeak_b) + np.array(ybshift )  

#  Convert peaks_piz to pixels
                peaks_pix = [[p1/0.2 for p1 in peak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in peak_b]]

                if v>0:
                    print " Arcsec: peaks_A = " , peak_a
                    print " Arcsec: peaks_B = " , peak_b
                    print " Pixels: peaks_pix = " ,  peaks_pix 

# Create the blended object using funct above 
#   blend = single image where we've added both
#   unblends = vector of 2 imgs, each of the ind ones
                blend, unblends = create_blend(peak_a, peak_b, e1a = e1ain,  e2a = e2ain, e1b = e1bin ,e2b = e2bin)

                if plotflag > presetval:
                    plt.title(" Img blended obj - (a+b) ")
                    if v>0:
                        print " >>>>>> Plotting blend.array - (unblends[0].array + unblends[1].array)  "
                    plt.imshow( blend.array - (unblends[0].array + unblends[1].array) , origin='lower');    
                    plt.colorbar()
                    plt.show()

# Use deblending code to separate them
#    templates = for each img
#    template_fractions
#    children = vector of 2 imgs, best estimates from deblending code
             #   templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix, interpolate=False, force_interpolate = False)
           #     templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix, interpolate=True, force_interpolate = False)  # Interpln on
                templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix, interpolate=True, force_interpolate = True)  # Force Interpln on
                if plotflag > presetval:
                    plt.title(" Template obj a ")
                    if v>0:
                        print " >>>>>> Plotting template fraction a"
                    plt.title(" Template 1")
                    plt.imshow( template_fractions[0] , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Template obj b ")
                    if v>0:
                        print " >>>>>> Plotting template fraction b"
                    plt.imshow( template_fractions[1] , origin='lower');    
                    plt.colorbar()
                    plt.show()

    ######### Plot children
                if plotflag > presetval:
                    if v>0:
                        print " >>>>>> Plotting Unblended img a "
                    plt.title(" Unblended img a ")
                    plt.imshow( unblends[0].array , origin='lower');            plt.colorbar();        plt.show()
                    if v>0:
                        print " >>>>>> Plotting Deblended child a "
                    plt.title(" Deblended child a ")
                    plt.imshow( children[0] , origin='lower');            plt.colorbar();        plt.show()
                    if v>0:
                        print " >>>>>> Plotting Resid of: (Deblended child a - Unblended img a) "
                    plt.title("Resid of: (Deblended child a - Unblended img a)  ")
                    plt.imshow( children[0] - unblends[0].array , origin='lower');            plt.colorbar();        plt.show()

                    if v>0:
                        print " >>>>>> Plotting Unblended img b "
                    plt.title(" Unblended img b ")
                    plt.imshow( unblends[1].array , origin='lower');            plt.colorbar();        plt.show()
                    if v>0:
                        print " >>>>>> Plotting Deblended child b "
                    plt.title(" Deblended child b ")
                    plt.imshow( children[1] , origin='lower');            plt.colorbar();        plt.show()
                    if v>0:
                        print " >>>>>> Plotting Resid of: (Deblended child b - Unblended img b) "
                    plt.title("Resid of: (Deblended child b - Unblended img b)  ")
                    plt.imshow( children[1] - unblends[1].array , origin='lower');            plt.colorbar();        plt.show()



    # Add the flux fraction imgs by hand 
    #   templatesum = adding two template_fraction imgs
                templatefracsum = template_fractions[0] + template_fractions[1]
    #    tarray = templatesum.array

                if plotflag > presetval:
                    plt.title(" Template fraction sum of a+b ")
                    if v>0:
                        print " >>>>>> Plotting templatesum"
                    plt.imshow( templatefracsum , origin='lower');    
                    plt.colorbar()
                    plt.show()


    # Add the child imgs by hand 
                sum_children = children[0] + children[1]    

                if plotflag > presetval:
                    plt.title(" Deblended children a + b")
                    if v>0:
                        print " >>>>>> Plotting sum_children "
                    plt.imshow( sum_children  , origin='lower');    
                    plt.colorbar()
                    plt.show()

        # Get resid
                residual = blend.array - sum_children

                if plotflag > presetval:
                    plt.title(" Residual: blended img - sum of deblended children ")
                    if v>0:
                        print " >>>>>> Plotting resid = blend.array - sum_children "
                    plt.imshow( residual  , origin='lower');    
                    plt.colorbar()
                    plt.show()

                if v>0:
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

                if v>0:
                    print "\n\n ********************* Fit to Unbl Obj a "
                    print "e1_a = ", e1_a, " ,  e1err = ", e1err
                    print "e2_a = ", e2_a, " ,  e2err = ", e2err

                e1a_unbl = e1_a ; e2a_unbl = e2_a 
                e1a_unblerr = e1err ; e2a_unblerr = e2err 

                # x and y posn
                x0_a = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_a = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit

                x0a_unbl = x0_a ;  y0a_unbl = y0_a

    ################# Unbl obj a plots
                if plotflag > (presetval-1):
                    plt.title(" Unblended fit results for a")
                    if v>0:
                        print " >>>>>> Plotting result of unbl fit obj a "
                    fitobj_a = drawgal(peak = peak_a , e1 = e1_a, e2 = e2_a)
                    plt.imshow( fitobj_a.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj a")
                    if v>0:
                        print " >>>>>> Plotting orig obj a "
                    plt.imshow( origimg.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (unblended fit result for a - orig obj a )")
                    if v>0:
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

                if v>0:
                    print "\n *********  Fit to Deblended Obj a "
                    print "e1_a = ", e1_a, " ,  e1err = ", e1err
                    print "e2_a = ", e2_a, " ,  e2err = ", e2err

                e1a_debl = e1_a ; e2a_debl = e2_a 
                e1a_deblerr = e1err ; e2a_deblerr = e2err 

                # x and y posn
                x0_a = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_a = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0a_debl = x0_a ;  y0a_debl = y0_a


    # Report the parameters to the interpreter screen                        
    #    lmfit.report_errors(mlresult.params)

    ################# Debl obj a plots
                if plotflag > (presetval-1):
                    plt.title(" Deblended fit results for a")
                    if v>0:
                        print " >>>>>> Plotting result of debl fit obj a "
                    fitobj_a = drawgal(peak = peak_a , e1 = e1_a, e2 = e2_a)
                    plt.imshow( fitobj_a.array  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Orig obj a")
                    if v>0:
                        print " >>>>>> Plotting orig obj a "
                    plt.imshow( origimg  , origin='lower');    
                    plt.colorbar()
                    plt.show()
                    plt.title(" Resid of:  (deblended fit result for a - orig obj a )")
                    if v>0:
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

                if v>0:
                    print "\n\n *********************  Fit to Unbl Obj b "
                    print "e1_b = ", e1_b, " ,  e1err = ", e1err
                    print "e2_b = ", e2_b, " ,  e2err = ", e2err
                    
                e1b_unbl = e1_b ; e2b_unbl = e2_b 
                e1b_unblerr = e1err ; e2b_unblerr = e2err 

                # x and y posn
                x0_b = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_b = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0b_unbl = x0_b ;  y0b_unbl = y0_b

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

                if v>0:
                    print "\n **********  Fit to Deblended Obj b "
                    print "e1_b = ", e1_b, " ,  e1err = ", e1err
                    print "e2_b = ", e2_b, " ,  e2err = ", e2err

                e1b_debl  = e1_b ; e2b_debl  = e2_b 
                e1b_deblerr = e1err ; e2b_deblerr = e2err 

                # x and y posn
                x0_b = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_b = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0b_debl = x0_b ;  y0b_debl = y0_b

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
                fitresults = [int(filenum), e1ain, e2ain, e1a_unbl,e1a_debl, e2a_unbl,e2a_debl,   e1bin, e2bin, e1b_unbl,e1b_debl, e2b_unbl, e2b_debl,   x0a_unbl,y0a_unbl, x0a_debl,y0a_debl, x0b_unbl,y0b_unbl, x0b_debl,y0b_debl ]
                fitresults.extend(peak_a)
                fitresults.extend(peak_b)

                fitdat.append(fitresults)
                
                if v>1:
                    print 'len(fitdat) = ', len(fitdat)

    ################################## End of all loops
    fitarray = np.array(fitdat)

    if v>1:
        print(fitarray)

    print " ----------------------------- sepstep = ", sepstep 

    if sepstep == 0:
        np.savetxt(outdir + 'deblendingTests_peak_A_'+str(origpeak_a) + '__peak_B_' + str(origpeak_b) +'_' + str(numfiles)+ '_runs.txt', fitarray, header="filenum   e1a_in e2a_in   e1a_unbl e1a_debl  e2a_unbl e2a_debl    e1b_in e2b_in  e1b_unbl  e1b_debl e2b_unbl e2b_debl    x0a_unbl y0a_unbl x0a_debl y0a_debl   x0b_unbl y0b_unbl   x0b_debl y0b_debl  x0_a y0_a  x0_b y0_b")
    else:
        np.savetxt(outdir + 'sepstep_' + str(sepstep/5), fitarray)
