
# MSSG, based on JEM and LCA code
# Start: 4/13/2015

# Three step proc:
# 1. simfit to A+B
# 2. deblend using output centers from the simfit
# 3. fit now to deblended objs

### Generic imports
import lmfit
#import ipdb
import sys
import random
# import triangle
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np

### Specific imports
import galsim
import mssg_deblend
import mssg_drawLibrary



################### Initze

#### Level to do printing at (setting it lower will print more stuff)
presetval = 1

initbasedev = galsim.BaseDeviate(1) # Random num seed -- when set to zero, uses machine time

########## Flags

### Flag to tell if if we want to do the simfit at all
simfit_flag = 0
 
#################################################################################### Function to plot gal figs
def makeplot(pltname, pltcontent, location = 'lower'):
          plt.title(pltname)
          print " >>>>>> Plotting ", pltname          
          dummyarr = np.array([1, 2, 3])
          if type(pltcontent) != type(dummyarr):
              print pltcontent, " is not a numpy array, will convert"
              pltcontent = pltcontent.array 
          plt.imshow( pltcontent , origin=location );    
          plt.colorbar()
          plt.show()
    
################################################################## Function to create and return simple single gal img
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
    print '\n In drawgal:  randnum = ', randnum

    img = convgal.drawImage(image=proto_image, method='phot', rng=randnum)
    img.array[np.where(img.array < 0)] = 0.
    return img




################################################################# Function to create the blended img
def create_blend(peak_a, peak_b, e1a = 0, e1b = 0 , e2a = 0, e2b = 0, rand = initbasedev, includepsf = False ):
    print '\n In create_blend:  randbasedev = ', rand
    # Create gaussian gal objs, sheared in various directions
    hlr_in = 1.0
    flux_in = 500 # 1e5 before
    gal1 = galsim.Gaussian(half_light_radius= hlr_in , flux= flux_in).shear(g1=e1a, g2= e2a).shift(peak_a)
    gal2 = galsim.Gaussian(half_light_radius= hlr_in , flux= flux_in).shear(g1=e1b, g2= e2b).shift(peak_b)
    
    # Add psf 
    psfshr = 0.00
    psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1 = psfshr,  e2 = -psfshr)

    convgal1 = galsim.Convolve([gal1,psf])
    convgal2 = galsim.Convolve([gal2,psf])
    
    # Now make imgs
    pixelsize = 49
    pixelscale = 0.2

    # Img1
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    image1 = convgal1.drawImage(image=proto_image, method='phot', rng=rand)
    #print '\n In drawgal:  randnum = ', rand
    image1.array[np.where(image1.array < 0)] = 0.
    if plotflag > presetval:
        plt.title(" FIRST PLOT:  Img obj a")
        plt.imshow( image1.array , origin='lower');    
        plt.colorbar()
        plt.show()
        print " >>>>>> Plotting img obj a"
        print " image1.array.sum() " , image1.array.sum()

    # Img2
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=rand)
    #print '\n In drawgal:  randnum = ', rand

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
    parser.add_argument("--prependname", default="", help="filename to prepend")
    parser.add_argument("--dirname", default="tmpdir/", help="output dir")
    parser.add_argument("--e1a", default=0, type=float, help="e1a in")
    parser.add_argument("--e2a", default=0, type=float, help="e2a in")
    parser.add_argument("--e1b", default=0, type=float, help="e1b in")
    parser.add_argument("--e2b", default=0, type=float, help="e2b in")
    parser.add_argument("--plotflag", default=0, type=int, help="Set to 1 to make plots")
    parser.add_argument("--centers", default='F', type=str, help="Set to E to use exact centers, F to use centers from simfit")
    parser.add_argument("--random", default=1, type=int, help="random number seed")

    args = parser.parse_args()

    outdir = args.dirname
    randseed = args.random # Getting the random number seed we'll pass on to create
    randbasedev=galsim.BaseDeviate(randseed) # Random num seed -- when set to zero, uses machine time, otw deterministic

    ########## Get ellips from cmd line
    e1ain = args.e1a
    e2ain = args.e2a

    e1bin = args.e1b
    e2bin = args.e2b

    plotflag = args.plotflag
    prename = args.prependname
    centers = args.centers  # Whether to use exact centers or to use centers from simfit when deblending

    # ************************************************************************ 4 arcsec sep
    # origpeak_a = (-2,0);   origpeak_b = (2,0)    

    # ************************************************************************ 2 arcsec sep
    origpeak_a = (-1,0);   origpeak_b = (1,0)    

    ################### Initze
    fitdat = []    # Vector that will have all output data in the end
    numfiles = 1  # Number of runs we want to do

    #### Normal input range i've been using
    e1a_range = [0.5,  0, -0.5]
    e1b_range = [0.5,  0, -0.5]

    #e1a_range = [0.5]
    #e1b_range = [0.5]

    ############################################################################## Begin loop
    for filenum in xrange(0,numfiles):

        ### Run over ellips
        for e1bin in e1b_range:
            for e1ain in e1a_range:

                print " \n\n ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin

                peak_a =  np.array(origpeak_a) ; peak_b =  np.array(origpeak_b) 
                print " \n\n\n peak_a = ",  peak_a 

                ########  Convert peaks_pix to pixels
                peaks_pix = [[p1/0.2 for p1 in peak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in peak_b]]

                print " Arcsec: peaks_A = " , peak_a
                print " Arcsec: peaks_B = " , peak_b
                print " Pixels: peaks_pix = " ,  peaks_pix 

                ####### Create the blended and unblended imgs    
                blend, unblends = create_blend(peak_a, peak_b, e1a = e1ain,  e2a = e2ain, e1b = e1bin ,e2b = e2bin,rand = randbasedev )

                if plotflag > presetval:
                    plt.title(" Img blended obj - (a+b) ")
                    print " >>>>>> Plotting blend.array - (unblends[0].array + unblends[1].array)  "
                    plt.imshow( blend.array - (unblends[0].array + unblends[1].array) , origin='lower');    
                    plt.colorbar()
                    plt.show()

                origpeak_a = peak_a ; origpeak_b = peak_b

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


                #################### Common vars to set
                galtype = galsim.Gaussian
                imsize = 49
                pixel_scale = 0.2
                dopsfconvln='y'
 
                if centers == 'F':
                
                    ############################################################################################### 1. Sim Fit
                    # This defines some seed values for when doing the lmfit object for galaxy one and two

                    # Parameters for object a
                    flux_a = 500 # 1e5 # 5e6          # total counts on the image
                    hlr_a = 1.0            # arcsec
                    e1_a = e1ain
                    e2_a = 0.0
                    x0_a = peak_a[0]
                    y0_a = peak_a[1]

                    # Parameters for object b
                    flux_b = flux_a       # total counts on the image
                    hlr_b = hlr_a         # arcsec
                    e1_b = e1bin
                    e2_b = 0.0
                    x0_b = peak_b[0]
                    y0_b = peak_b[1]

                    # Declare an array of all the seed values
                    p0 = 1.0*np.array([flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                                       flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b])  # These are just the init guesses for the fit (all pre-defined nums from above)

                    print " \n ************************************ 1. About to do first simfit"


                    ############################################## Sim Fit to Obj a+b
                    params = lmfit.Parameters()
                    params.add('flux_a', value=p0[0], min=0.0)   
                    params.add('hlr_a', value=p0[1], min=0.0)
                    params.add('e1_a', value=p0[2], min=-1.0, max=1.0)
                    params.add('e2_a', value=p0[3], min=-1.0, max=1.0)
                    params.add('x0_a',value=p0[4])
                    params.add('y0_a',value=p0[5])

                    params.add('flux_b', value=p0[6], min=0.0)
                    params.add('hlr_b', value=p0[7], min=0.0)
                    params.add('e1_b', value=p0[8], min=-1.0, max=1.0)
                    params.add('e2_b', value=p0[9], min=-1.0, max=1.0)
                    params.add('x0_b',value=p0[10])
                    params.add('y0_b',value=p0[11])

                    tot_unbl_img =  unblends[0]+unblends[1]
                    #ipdb.set_trace()

                    result = lmfit.minimize(mssg_drawLibrary.resid_2obj,   params,   args=( tot_unbl_img , imsize,imsize,pixel_scale, galtype, galtype ))

                    # Report the parameters to the interpreter screen                        
                    #lmfit.report_errors(result.params)

                    #################### Extract the centers from the simfit
                    x0_a_sf = result.params['x0_a'].value
                    y0_a_sf = result.params['y0_a'].value

                    x0_b_sf = result.params['x0_b'].value
                    y0_b_sf = result.params['y0_b'].value

                    e1_a_sf = result.params['e1_a'].value  # Get out e1 val of obj a from fit
                    e1_b_sf = result.params['e1_b'].value  # Get out e1 val of obj b from fit
                
                    print "\n\n\n ************** In orig simfit:  x0_a_sf ,  y0_a_sf = ", x0_a_sf ,  y0_a_sf
                    print " x0_b_sf ,  y0_b_sf = ", x0_b_sf ,  y0_b_sf


                    #    ipdb.set_trace()


            ##################################################################################################### 2. Deblend
            # Use deblending code to separate them
            #    templates = for each img
            #    template_fractions
            #    children = vector of 2 imgs, best estimates from deblending code

                print " \n ***************************************** 2. About to deblend"

               ########################### Which centers to use
                if centers == 'F':
               ########################### Using fit vals as centers
                    print 'Using SimFit centers'
                    curpeak_a = (x0_a_sf, y0_a_sf);       curpeak_b = (x0_b_sf, y0_b_sf)       
                else:
               ########################### Using known exact vals as centers
                    print 'Using Exact centers'
                    curpeak_a = origpeak_a;       curpeak_b = origpeak_b
                
#                print "\n\n\n x0_a_sf ,  y0_a_sf = ", x0_a_sf ,  y0_a_sf
                print "    curpeak_a,  curpeak_b = ",  curpeak_a,  curpeak_b                              
                #  Convert peaks_pix to pixels
                peaks_pix = [[p1/0.2 for p1 in curpeak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in curpeak_b]]

                origpeaks_pix = [[p1/0.2 for p1 in origpeak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in origpeak_b]]
                
                print " peaks_pix ,  origpeaks_pix = ", peaks_pix ,  origpeaks_pix 

                templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix, interpolate=False, force_interpolate = False)

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                # origtemplates, origtemplate_fractions, origchildren = mssg_deblend.deblend(blend.array, origpeaks_pix, interpolate=False, force_interpolate = False)
#                ipdb.set_trace()


                ########## Plot template
                if plotflag > presetval:
                    pltname = " Template obj a "
                    pltcontent = template_fractions[0] 
                    makeplot(pltname, pltcontent)

                    ######### Plot children
                    ### Obj a
                    pltname = " Unblended img a ";        pltcontent = unblends[0]
                    makeplot(pltname, pltcontent)

                    pltname = " Deblended img a ";        pltcontent = children[0]
                    makeplot(pltname, pltcontent)

                    pltname = " Resid of: (Deblended child a - Unblended img a) ";        pltcontent = children[0] - unblends[0].array 
                    makeplot(pltname, pltcontent)

            ### Obj b
                    pltname = " Unblended img b ";        pltcontent = unblends[1]
                    makeplot(pltname, pltcontent)

                    pltname = " Deblended img b ";        pltcontent = children[1]
                    makeplot(pltname, pltcontent)

                    pltname = " Resid of: (Deblended child b - Unblended img b) ";        pltcontent = children[1] - unblends[1].array 
                    makeplot(pltname, pltcontent)






                ################################################################## 3. Run deblended fits
                print " \n ************************************** 3. Now about to run lmfit for final separate deblended fits "
                # Common params to all



               ########################### Which centers to use
                if centers == 'F':
               ########################### Using fit vals as centers
                    print 'Using SimFit centers'
                    curpeak_a = (x0_a_sf, y0_a_sf);       curpeak_b = (x0_b_sf, y0_b_sf)       
                else:
               ########################### Using known exact vals as centers
                    print 'Using Exact centers'
                    curpeak_a = origpeak_a;       curpeak_b = origpeak_a


                peak_a =  np.array(curpeak_a) ; peak_b =  np.array(curpeak_b) 
                #    print " \n\n\n peak_a = ",  peak_a 

                #  Convert peaks_pix to pixels
                peaks_pix = [[p1/0.2 for p1 in peak_a],  # Div by 0.2 to convert back to pixels
                                         [p2/0.2 for p2 in peak_b]]

                print " \n\n\n  ************************* Post deblending "   
                print " Arcsec: peaks_A = " , peak_a
                print " Arcsec: peaks_B = " , peak_b
                print " Pixels: peaks_pix = " ,  peaks_pix 

               # The below are just initial guesses for the fitter
                minmaxval = 0.707 # Temporary limit for the fitter till i fix it better by doing value in quadrature
                fit_params = lmfit.Parameters()  
                fit_params.add('e1_a', value=e1ain, min= -minmaxval, max=minmaxval)
                fit_params.add('e2_a', value=0.0, min= -minmaxval, max=minmaxval)
                fit_params.add('flux_a', value=2000.0)
                fit_params.add('x0_a', value=peak_a[0])
                fit_params.add('y0_a', value=peak_a[1])
                fit_params.add('hlr_a', value=1)

            ####### Whether to convolve with PSF
                dopsfconvln = 'y'

                ###################################### Obj a
                print " ************************* Beginning fits"
                 
                 ############# Unbl Obj a        
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

                        # x and y posn
                x0_a = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_a = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit

                x0a_unbl = x0_a ;  y0a_unbl = y0_a


                ############ Deblended Obj a
                origimg = children[0]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) ) 

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #                origkid0 = origchildren[0]    
                #                origmlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origkid0,imsize,imsize,pixel_scale, galtype, dopsfconvln) ) 

                # Report the parameters to the interpreter screen                        
                print "********************* Post second fit -- obj A"
                #lmfit.report_errors(mlresult.params)

                # Extract vals
                e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #orige1_a = origmlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                #orige2_a = origmlresult.params['e2_a'].value  # Get out e2 val of obj a from fit


                print "\n *********  Deblended Obj a "
                print "e1_a = ", e1_a, " ,  e1err = ", e1err
                print "e2_a = ", e2_a, " ,  e2err = ", e2err

                e1a_debl = e1_a ; e2a_debl = e2_a 
                e1a_deblerr = e1err ; e2a_deblerr = e2err 

                # x and y posn
                x0_a = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_a = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0a_debl = x0_a ;  y0a_debl = y0_a


                ########### Unbl Obj b        
                fit_params.add('e1_a', value=e1bin, min= -minmaxval, max=minmaxval)  # Set to actual input value for obj b
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

                # x and y posn
                x0_b = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_b = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0b_unbl = x0_b ;  y0b_unbl = y0_b

                ######## Deblended Obj b
                origimg = children[1]    
                mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) )  

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #origkid1 = origchildren[1]    
                #origmlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origkid1,imsize,imsize,pixel_scale, galtype, dopsfconvln) ) 


                # Report the parameters to the interpreter screen                        
                print "********************* Post second fit -- obj B"
                #lmfit.report_errors(mlresult.params)

               # Extract vals  -- note that calling below: mlresult.params['e1_a']  is correct because the '_a' name is defined that way inside the fitter
                e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
                e1err = np.sqrt(np.diag(mlresult.covar)[0])
                e2err = np.sqrt(np.diag(mlresult.covar)[1])

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #orige1_b = origmlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
                #orige2_b = origmlresult.params['e2_a'].value  # Get out e2 val of obj a from fit

                print "\n **********  Deblended Obj b "
                print "e1_b = ", e1_b, " ,  e1err = ", e1err
                print "e2_b = ", e2_b, " ,  e2err = ", e2err

                e1b_debl  = e1_b ; e2b_debl  = e2_b 
                e1b_deblerr = e1err ; e2b_deblerr = e2err 

                # x and y posn
                x0_b = mlresult.params['x0_a'].value  # Get out x0 val of obj a from fit
                y0_b = mlresult.params['y0_a'].value  # Get out y0 val of obj a from fit
                x0b_debl = x0_b ;  y0b_debl = y0_b


                ########################## Final results
                print "\n\n\n \n **********  Final --- Deblended Obj a "
                print "e1_a = ", e1_a, " e1ain - e1_a = ", e1ain - e1_a 

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #print "orige1_a = ", e1_a, " e1ain - orige1_a = ", e1ain - orige1_a 

                print "\n **********  Deblended Obj b "
                print "e1_b = ", e1_b, " e1bin - e1_b = ", e1bin - e1_b 

                #### If testing to compare origtemplates -- i.e. those with exactcenters
                #print "orige1_b = ", e1_b, " e1bin - orige1_b = ", e1bin - orige1_b 

                    ################### Result vec for this fit

                # Case of simfit centers
                if centers == 'F':
                    fitresults = [int(filenum), e1ain, e2ain, e1a_unbl,e1a_debl, e2a_unbl,e2a_debl,   e1bin, e2bin, e1b_unbl,e1b_debl, e2b_unbl, e2b_debl,   x0a_unbl,y0a_unbl, x0a_debl,y0a_debl, x0b_unbl,y0b_unbl, x0b_debl,y0b_debl ,x0_a_sf, x0_b_sf, e1_a_sf, e1_b_sf]

                # In case of exact centers, we don't want to refer to any of the simfit outputs
                fitresults = [int(filenum), e1ain, e2ain, e1a_unbl,e1a_debl, e2a_unbl,e2a_debl,   e1bin, e2bin, e1b_unbl,e1b_debl, e2b_unbl, e2b_debl,   x0a_unbl,y0a_unbl, x0a_debl,y0a_debl, x0b_unbl,y0b_unbl, x0b_debl,y0b_debl ] 

                fitresults.extend(peak_a)
                fitresults.extend(peak_b)

                fitdat.append(fitresults)

                print 'len(fitdat) = ', len(fitdat)

    ################################## End of all loops
    fitarray = np.array(fitdat)

    #print(fitarray)
    
    #np.savetxt(outdir+prename+'deblendingTests_peak_A_'+str(origpeak_a) + '__peak_B_' + str(origpeak_b) +  '_runnum_'+ str(randseed) +'.txt', fitarray, header="filenum   e1a_in e2a_in   e1a_unbl e1a_debl  e2a_unbl e2a_debl    e1b_in e2b_in  e1b_unbl  e1b_debl e2b_unbl e2b_debl    x0a_unbl y0a_unbl x0a_debl y0a_debl   x0b_unbl y0b_unbl   x0b_debl y0b_debl  x0_a y0_a  x0_b y0_b")

    np.savetxt(outdir+prename+'tmpfile_runnum_'+ str(randseed) +'.txt', fitarray)
