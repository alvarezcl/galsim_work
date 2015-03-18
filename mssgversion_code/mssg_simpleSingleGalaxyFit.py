# Modifications by MSSG
# Start: 3/17/2015
# based on:  mssg_2D_Sersic_vary_distance.py

########################## Make simple Gaussian galaxy and fit it -- to get the technique down

#### Generic
from __future__ import division
import lmfit
import numpy as np
import ipdb
import sys
import matplotlib.gridspec as gridspec

#### Specific
import galsim
import mssg_drawLibrary
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- Begin with parameters for the image 

# Parameters for object a
flux_a = 5e6          # total counts on the image
HLR_a = 3.            # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = 0  # We'll step gal A to the left in steps of 0.5*step_sep - mg (Shift in arcsec from center)
y0_a = 0


# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 0.2     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
imsize = 5*HLR_a/pixel_scale            # Total size of image on a side, in pixels - i changed name from size to imsize -mg
dopsfconvln = 'y'
x_len = y_len = imsize




# Galsim function definitions 
galtype = galsim.Gaussian

# Define a numpy array of params - mg
param_array = np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])

    
# ------------------------------------------------------------------------ Create the image 
    
# Random Seed
# Putting in a zero in the arg means it will vary from plot to plot, *not* exactly same each time - mg
randomnum = galsim.BaseDeviate(0) 
    
# Fill in the images - mg

img = mssg_drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                         imsize,imsize,pixel_scale,galtype,randomnum)
     

# Obtain the image bounds and domain information
x_cen,y_cen,x,y,X,Y = mssg_drawLibrary.return_domain(img)
H = img.array # I believe this returns the 2D numpy array assoc with the image - mg

p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])  # These are all pre-defined nums from above - mg

parameters = lmfit.Parameters()

parameters.add('flux_a', value=p0[0])   
parameters.add('hlr_a', value=p0[1], min=0.0)
parameters.add('e1_a', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2_a', value=p0[3], min=-1.0, max=1.0)
parameters.add('x0_a',value=p0[4])
parameters.add('y0_a',value=p0[5])

result = lmfit.minimize(mssg_drawLibrary.resid_1obj,   parameters,   args=(img ,imsize,imsize,pixel_scale, galtype, galtype ))
        
lmfit.report_errors(result.params)


best_fit = mssg_drawLibrary.draw_galaxy_1(result.params['flux_a'].value,
                                           result.params['hlr_a'].value,
                                           result.params['e1_a'].value,
                                           result.params['e2_a'].value,
                                           result.params['x0_a'].value,
                                           result.params['y0_a'].value, x_len,y_len, pixel_scale, galtype, dopsfconvln)


############################################ Set pltflag
pltflg = 1

if pltflg == 1:
    # Initialize the geometry of the grid
    gs = gridspec.GridSpec(1,3)
    fig = plt.figure(figsize=(10,5))
    # Plot the orig data
    f1 = fig.add_subplot(gs[0,0])
    plt.title('Binned Image of Galaxy')
    fig1 = f1.imshow(H,interpolation='none',origin='lower')
    plt.colorbar(fig1, shrink=0.5) 

    # Plot the best fit 
    f2 = fig.add_subplot(gs[0,1])
    plt.title('Best Fit')
    fig2 = f2.imshow(best_fit.array,interpolation='none')
    plt.colorbar(fig2, shrink=0.5)

    # Plot the Resid
    f3 = fig.add_subplot(gs[0,2])
    plt.title('Resid')
    fig3 = f3.imshow(best_fit.array - H,interpolation='none')
    plt.colorbar(fig3, shrink=0.5)
    plt.show()
    
sys.exit()


