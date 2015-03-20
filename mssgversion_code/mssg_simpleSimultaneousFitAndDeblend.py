# MSSG, based on JEM and LCA code
# Start: 3/10/2015

### Generic imports
import lmfit
import ipdb
import sys
import random
import triangle
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np

### Specific imports
import galsim
import mssg_deblend
import mssg_drawLibrary



################### Initze

#### Level to do printing at (setting it lower will print more stuff)
presetval = 0

randnum=galsim.BaseDeviate(0) # Random num seed -- when set to zero, uses machine time



#################################################### Function to draw and return simple single gal img
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
    img = convgal.drawImage(image=proto_image, method='phot', rng=randnum)
    img.array[np.where(img.array < 0)] = 0.
    return img




###################################################### Function to create the blended img
def create_blend(peak_a, peak_b, e1a = 0, e1b = 0 , e2a = 0, e2b = 0):
    # Create gaussian gal objs, sheared in various directions
    gal1 = galsim.Gaussian(fwhm=1.0, flux=100000.0).shear(g1=e1a, g2= e2a).shift(peak_a)
    gal2 = galsim.Gaussian(fwhm=1.0, flux=100000.0).shear(g1=e1b, g2= e2b).shift(peak_b)
    
    # Add psf 
    psfshr = 0.00
    psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1=psfshr,  e2=-0.0)

    convgal1 = galsim.Convolve([gal1,psf])
    convgal2 = galsim.Convolve([gal2,psf])
    
    # Now make imgs
    pixelsize = 49
    pixelscale = 0.2

    # Img1
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    image1 = convgal1.drawImage(image=proto_image, method='phot', rng=randnum)
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
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=randnum)
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
    e1ain = args.e1a
    e2ain = args.e2a

    e1bin = args.e1b
    e2bin = args.e2b

    plotflag = args.plotflag
    origpeak_a = (-2,0);   origpeak_b = (2,0)    # Horiz sep - centers separated by 1.6", to match with Luis, EXACTLY

    fitdat = []

    peak_a =  np.array(origpeak_a) ; peak_b =  np.array(origpeak_b) 
    print " \n\n\n peak_a = ",  peak_a 

#  Convert peaks_pix to pixels
    peaks_pix = [[p1/0.2 for p1 in peak_a],  # Div by 0.2 to convert back to pixels
                             [p2/0.2 for p2 in peak_b]]
    
    print " Arcsec: peaks_A = " , peak_a
    print " Arcsec: peaks_B = " , peak_b
    print " Pixels: peaks_pix = " ,  peaks_pix 
    
    blend, unblends = create_blend(peak_a, peak_b, e1a = e1ain,  e2a = e2ain, e1b = e1bin ,e2b = e2bin)

    if plotflag > presetval:
        plt.title(" Img blended obj - (a+b) ")
        print " >>>>>> Plotting blend.array - (unblends[0].array + unblends[1].array)  "
        plt.imshow( blend.array - (unblends[0].array + unblends[1].array) , origin='lower');    
        plt.colorbar()
        plt.show()

######################################################### Sim Fit
# Parameters for object a
    flux_a = 5e6          # total counts on the image
    HLR_a = 3.            # arcsec
    e1_a = 0.0
    e2_a = 0.0
    x0_a = peak_a[0]
    y0_a = peak_a[1]

# Parameters for object b
    flux_b = flux_a       # total counts on the image
    HLR_b = HLR_a         # arcsec
    e1_b = 0.0
    e2_b = 0.0
    x0_b = peak_b[0]
    y0_b = peak_b[1]

# Define some seed that's far from true values and insert into
# lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                       flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b]) # These are all pre-defined nums from above - mg
    params = lmfit.Parameters()
# Shouldn't these really be all called _a and then _b vs. 1 and 2..?  Maybe they're just name labels though..?  -mg
    params.add('flux_a', value=p0[0])   
    params.add('hlr_a', value=p0[1], min=0.0)
    params.add('e1_a', value=p0[2], min=-1.0, max=1.0)
    params.add('e2_a', value=p0[3], min=-1.0, max=1.0)
    params.add('x0_a',value=p0[4])
    params.add('y0_a',value=p0[5])

    params.add('flux_b', value=p0[6])
    params.add('hlr_b', value=p0[7], min=0.0)
    params.add('e1_b', value=p0[8], min=-1.0, max=1.0)
    params.add('e2_b', value=p0[9], min=-1.0, max=1.0)
    params.add('x0_b',value=p0[4])
    params.add('y0_b',value=p0[5])

    galtype = galsim.Gaussian
    imsize = 49
    pixel_scale = 0.2


    print " ************** About to fit"

    result = lmfit.minimize(mssg_drawLibrary.resid_2obj,   params,   args=(blend, imsize,imsize,pixel_scale, galtype, galtype ))

    ipdb.set_trace()

# Report the parameters to the interpreter screen                        
    lmfit.report_errors(result.params)

         
######################################################### Deblend
# Use deblending code to separate them
#    templates = for each img
#    template_fractions
#    children = vector of 2 imgs, best estimates from deblending code

# ....................
