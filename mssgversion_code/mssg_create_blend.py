import matplotlib.pyplot as plt
import numpy as np
import galsim
import random

###################################################### Function to create the blended img
def create_blend(peak_a, peak_b, e1a = 0, e1b = 0 , e2a = 0, e2b = 0, imgsize = 0, pixelscale = 0, mrkrsize = 10):

    # Global vars
    randnum=galsim.BaseDeviate(1) # Random num seed -- when set to zero, uses machine time

    #### Level to do printing at (setting it lower will print more stuff)
    presetval = 0
    #    plotflag =  0  # The default, to plot the minimum
    plotflag =  1  # Print more
    
    # Create gaussian gal objs, sheared in various directions
    hlr_in = 1.0
    flux_in = 1e5
    gal1 = galsim.Gaussian(half_light_radius= hlr_in , flux= flux_in).shear(g1=e1a, g2= e2a).shift(peak_a)
    gal2 = galsim.Gaussian(half_light_radius= hlr_in , flux= flux_in/1.2 ).shear(g1=e1b, g2= e2b).shift(peak_b)
    
    # Add psf 
    psfshr = 0.00
    psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1 = psfshr,  e2 = -psfshr)

    convgal1 = galsim.Convolve([gal1,psf])
    convgal2 = galsim.Convolve([gal2,psf])
    
    # Now make imgs
#    imgsize = 49
#    pixelscale = 0.2
    imgcent = imgsize/2

    print "peak_a[0] = ", peak_a[0]
    print 'imgcent = ', imgcent
    acent = (imgcent+peak_a[0]/pixelscale, imgcent+peak_a[1]/pixelscale)
    bcent = (imgcent+peak_b[0]/pixelscale, imgcent+peak_b[1]/pixelscale)

    # Img1
    proto_image = galsim.ImageD(imgsize, imgsize, scale = pixelscale)
    image1 = convgal1.drawImage(image=proto_image, method='phot', rng=randnum)
#    image1 = convgal1.drawImage(image=proto_image)
    print "************** About to print where img pixels are < 0"
    print image1.array[np.where(image1.array < 0)] 
    image1.array[np.where(image1.array < 0)] = 0.
    if plotflag > presetval:
        plt.title(" FIRST PLOT:  Img obj a")
        plt.imshow( image1.array , origin='lower', interpolation='none');    
        plt.colorbar()
        plt.plot( acent[0], acent[1], color='g', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
        plt.show()
        print " >>>>>> Plotting img obj a"
        print " image1.array.sum() " , image1.array.sum()

    # Img2
    proto_image = galsim.ImageD(imgsize, imgsize, scale = pixelscale)
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=randnum)
#    image2 = convgal2.drawImage(image=proto_image)
    image2.array[np.where(image2.array < 0)] = 0.    
    if plotflag > presetval: 
        plt.title(" Img obj b")
        print " >>>>>> Plotting img obj b"
        plt.imshow( image2.array , origin='lower',  interpolation='none' );    
        plt.colorbar()
        plt.plot(bcent[0], bcent[1],  color='y', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
        plt.show()

### Add them into one image
    imagesum =  image1+image2

    if plotflag > presetval:
        plt.title(" Img obj a + b")
        print " >>>>>> Plotting imgsum"
        plt.plot( acent[0], acent[1], color='g', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
        plt.plot(bcent[0], bcent[1],  color='y', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
        
        plt.imshow( imagesum.array , origin='lower',  interpolation='none' );    
        plt.colorbar()
        plt.show()

        
    return imagesum, [image1, image2]




