# MSSG, based on JEM and LCA code
# 10/15/2014

import galsim
import numpy as np
import mssg_deblend
import mssg_drawLibrary
import lmfit
import ipdb
from argparse import ArgumentParser
import matplotlib.pyplot as plt

bd=galsim.BaseDeviate(2) # Random num seed

############## Function to create the blended img

def create_blend(peak1, peak2, e1a = 0, e1b = 0 , e2a = 0, e2b = 0):
    # Create gaussian gal objs, sheared in various directions


    gal1 = galsim.Gaussian(fwhm=1.0, flux=2500.0).shear(g1=e1a, g2= e2a).shift(peak1)
    gal2 = galsim.Gaussian(fwhm=1.0, flux=2500.0).shear(g1=e1b, g2= e2b).shift(peak2)

    
    # Add psf 
    psfshr = 0.00
    psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1=psfshr,  e2=-0.0)

    convgal1 = galsim.Convolve([gal1,psf])
    convgal2 = galsim.Convolve([gal2,psf])
    
    convgal1 = gal1 
    convgal2 = gal2


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

 #    plt.imshow( image1.array , origin='lower');    
#    plt.show()

    # Img2
    proto_image = galsim.ImageD(pixelsize, pixelsize, scale = pixelscale)
    #    image2 = gal2.drawImage(image=proto_image, method='fft')
    #image2 = gal2.drawImage(image=proto_image, method='phot', rng=bd)
    # With psf convln
    #image2 = convgal2.drawImage(image=proto_image, method='fft')
    image2 = convgal2.drawImage(image=proto_image, method='phot', rng=bd)
    image2.array[np.where(image2.array < 0)] = 0.
    

#    plt.imshow( image2.array , origin='lower');    
#    plt.show()

#    ipdb.set_trace()
    return image1+image2, [image1, image2]




############################################################################# Main
if __name__ == '__main__':
# Parse command line args
    parser = ArgumentParser()
    parser.add_argument("--outfile", default="deblendsOutput/deblendingTests", help="output text filename")
    parser.add_argument("--e1a", default=0, type=float, help="e1a in")
    parser.add_argument("--e2a", default=0, type=float, help="e2a in")
    parser.add_argument("--e1b", default=0, type=float, help="e1b in")
    parser.add_argument("--e2b", default=0, type=float, help="e2b in")

    args = parser.parse_args()

# Get ellips from cmd line
    e1a_in = args.e1a
    e2a_in = args.e2a

    e1b_in = args.e1b
    e2b_in = args.e2b

# Centroids
    peak1 = (-1, 0)
    peak2 = (+1, 0)

    peaks_pix = [[p1/0.2 for p1 in peak1],
                 [p2/0.2 for p2 in peak2]]


    blend, unblends = create_blend(peak1, peak2, e1a = e1a_in,  e2a = e2a_in, e1b = e1b_in ,e2b = e2b_in)
    templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix)

    templatesum = template_fractions[0] + template_fractions[1]
#    tarray = templatesum.array

    sum_children = children[0] + children[1]    
    residual = blend.array - sum_children

    print 'blend.array.sum() = ' , blend.array.sum()
    print 'children[0].sum() + children[1].sum() = ', children[0].sum() + children[1].sum()

    # Run fits
    # Common params
    print " ***** Now about to run lmfit "
    imsize = 49
    pixel_scale = 0.2
    galtype = galsim.Gaussian

    minmaxval = 0.7 # Temporary till i fix it better by doing value in quadrature
    fit_params = lmfit.Parameters()  # The below are just initial guesses for the fitter
    fit_params.add('e1_a', value=0.1, min= -minmaxval, max=minmaxval)
    fit_params.add('e2_a', value=0.3, min= -minmaxval, max=minmaxval)
    fit_params.add('flux_a', value=2000.0)
    fit_params.add('x0_a', value=0)
    fit_params.add('y0_a', value=0)
    fit_params.add('hlr_a', value=0.43)

## Exact Obj a        
    origimg = unblends[0]    
    dopsfconvln = 'n'
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype , dopsfconvln) ) 
    
    # Extract vals
    e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
    
    print "\n\n ********************* Exact Obj a "
    print "e1_a = ", e1_a
    print "e2_a = ", e2_a
    
    e1a_exact = e1_a ; e2a_exact = e2_a 

    # Report the parameters to the interpreter screen                        
    #    lmfit.report_errors(mlresult.params)

## Deblended Obj a
    origimg = children[0]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) ) 
    
# Extract vals
    e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
    
    print "\n *********  Deblended Obj a "
    print "e1_a = ", e1_a
    print "e2_a = ", e2_a
    
    e1a_debl = e1_a ; e2a_debl = e2_a 

# Report the parameters to the interpreter screen                        
#    lmfit.report_errors(mlresult.params)

## Obj b        
    origimg = unblends[1]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) )  
    
    # Extract vals
    e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
    
    print "\n\n *********************  Exact Obj b "
    print "e1_b = ", e1_b
    print "e2_b = ", e2_b
    
    e1b_exact = e1_b ; e2b_exact = e2_b 

    # Report the parameters to the interpreter screen                        
    #    lmfit.report_errors(mlresult.params)

## Deblended Obj b
    origimg = children[1]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype, dopsfconvln) )  
    
    # Extract vals  -- note that calling below: mlresult.params['e1_a']  is correct because the '_a' name is defined that way inside the fitter
    e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
        
    print "\n **********  Deblended Obj b "
    print "e1_b = ", e1_b
    print "e2_b = ", e2_b
    
    e1b_debl  = e1_b ; e2b_debl  = e2_b 

    # Report the parameters to the interpreter screen                        
#    lmfit.report_errors(mlresult.params)
              
        #            HLRvec.append(ml.params['hlr'].value)
        #            fluxvec.append(ml.params['flux'].value)
    

# open(args.outfile+"_e1aIn_"+str(e1a_in)+"_e2aIn_"+str(e2a_in)+".txt", 'wa') as f:

    with open(args.outfile+".txt", 'a') as f:
        f.write("{}\t".format(e1a_in))
        f.write("{}\t".format(e1a_exact))
        f.write("{}\t".format(e1a_debl) )
        f.write("{}\n".format(e1a_debl - e1a_exact) )
#        f.write("{}\t".format(e2a_in), "{}\t".format(e2a_exact), "{}\n".format(e2a_debl) )


'''
################################################## Plot figs
    fig = plt.figure(figsize=(20,11))
    ax1 = fig.add_subplot(4,3,1)

    # Orig imgs
    plt.title('Original Image')
    a = ax1.imshow(blend.array, vmin=0, vmax=50)
    plt.colorbar(a)
    ax2 = fig.add_subplot(4,3,2)
    plt.title('Object 1')
    ax3 = fig.add_subplot(4,3,3)
    plt.title('Object 2')  
    ax2.imshow(unblends[0].array, vmin=0, vmax=50)
    ax3.imshow(unblends[1].array, vmin=0, vmax=50)

    # Symmetrized 'templates'
    ax5 = fig.add_subplot(4,3,5)
    plt.title('Symmetric Template about Peak 1')
    ax6 = fig.add_subplot(4,3,6)
    plt.title('Symmetric Template about Peak 2')
    ax5.imshow(templates[0], vmin=0, vmax=50)
    ax6.imshow(templates[1], vmin=0, vmax=50)

    # Template fractions
    ax8 = fig.add_subplot(4,3,8)
    plt.title('Template Fraction of Object 1')
    ax9 = fig.add_subplot(4,3,9)
    plt.title('Template Fraction of Object 2')    
    ax8.imshow(template_fractions[0], vmin=0, vmax=1)
    ax9.imshow(template_fractions[1], vmin=0, vmax=1)
    
    # Individual children
    ax11 = fig.add_subplot(4,3,11)
    plt.title('Deblended Object 1')
    ax12 = fig.add_subplot(4,3,12)
    plt.title('Deblended Object 2')
    ax11.imshow(children[0], vmin=0, vmax=50)
    ax12.imshow(children[1], vmin=0, vmax=50)

    # Resid    
    fig = plt.figure(2)
    ax1 = fig.add_subplot(1,1,1)
    a = ax1.imshow(residual,interpolation='none')
    plt.colorbar(a)

    # Resid    
    fig = plt.figure(3)
    ax1 = fig.add_subplot(1,1,1)
    a = ax1.imshow(templatesum,interpolation='none')
    plt.colorbar(a)

    plt.show()
'''



