# MSSG, based on JEM and LCA code
# 10/15/2014

import galsim
import numpy as np
import mssg_deblend
import mssg_drawLibrary
import lmfit
import ipdb

bd=galsim.BaseDeviate(1) # Random num seed

def create_blend(peak1, peak2):
    # Create gaussian gal objs, sheared in various directions
    gal1 = galsim.Gaussian(fwhm=1.2, flux=2000.0).shear(g1=0.1, g2=0.3).shift(peak1)
    gal2 = galsim.Gaussian(fwhm=1.8, flux=2500.0).shear(g1=-0.1, g2=-0.4).shift(peak2)

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

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    peak1 = (0.8, 0.2)
    peak2 = (-0.3, 1.26)

    peaks_pix = [[p1/0.2 for p1 in peak1],
                 [p2/0.2 for p2 in peak2]]

    blend, unblends = create_blend(peak1, peak2)
    templates, template_fractions, children = mssg_deblend.deblend(blend.array, peaks_pix)
    sum_children = children[0] + children[1]    
    residual = blend.array - sum_children

    print blend.array.sum()
    print children[0].sum() + children[1].sum()

    # Run fits
    # Common params
    print " ***** Now about to run lmfit "
    imsize = 49
    pixel_scale = 0.2
    galtype = galsim.Gaussian

    minmaxval = 0.7 # Temporary till i fix it better by doing value in quadrature
    fit_params = lmfit.Parameters()  # The below are just initial guesses for the fitter
    fit_params.add('e1_a', value=0.3, min= -minmaxval, max=minmaxval)
    fit_params.add('e2_a', value=0.2, min= -minmaxval, max=minmaxval)
    fit_params.add('flux_a', value=2000.0)
    fit_params.add('x0_a', value=0)
    fit_params.add('y0_a', value=0)
    fit_params.add('hlr_a', value=0.43)

## Exact Obj a        
    origimg = unblends[0]    
    dopsfconvln = 'y'
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype , dopsfconvln) ) 
    
        # Extract vals
    e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
        
    print "\n\n ********************* Exact Obj a "
    print "e1_a = ", e1_a
    print "e2_a = ", e2_a
    
    # Report the parameters to the interpreter screen                        
    lmfit.report_errors(mlresult.params)

## Deblended Obj a
    origimg = children[0]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype ) ) 
    
# Extract vals
    e1_a = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_a = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
        
    print "\n\n *********************  Deblended Obj a "
    print "e1_a = ", e1_a
    print "e2_a = ", e2_a
    
# Report the parameters to the interpreter screen                        
    lmfit.report_errors(mlresult.params)

## Obj b        
    origimg = unblends[1]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype ) ) 
    
        # Extract vals
    e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
        
    print "\n\n *********************  Exact Obj b "
    print "e1_b = ", e1_b
    print "e2_b = ", e2_b
    
    # Report the parameters to the interpreter screen                        
    lmfit.report_errors(mlresult.params)

## Deblended Obj 1
    origimg = children[1]    
    mlresult = lmfit.minimize(mssg_drawLibrary.resid_1obj, fit_params, args=(origimg,imsize,imsize,pixel_scale, galtype ) ) 
    
        # Extract vals
    e1_b = mlresult.params['e1_a'].value  # Get out e1 val of obj a from fit
    e2_b = mlresult.params['e2_a'].value  # Get out e2 val of obj a from fit
        
    print "\n\n *********************  Deblended Obj b "
    print "e1_b = ", e1_b
    print "e2_b = ", e2_b
    
    # Report the parameters to the interpreter screen                        
    lmfit.report_errors(mlresult.params)

        
        
        #            HLRvec.append(ml.params['hlr'].value)
        #            fluxvec.append(ml.params['flux'].value)
    


# '''
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

    plt.show()
# '''
