# Modifications by MSSG
# Start: 3/17/2015
# based on:  mssg_2D_Sersic_vary_distance.py

################################################################# Make simple Gaussian galaxy and fit it -- fully commented

#### Generic imports
import lmfit
import numpy as np
import ipdb
import sys
import matplotlib.gridspec as gridspec

#### Specific imports
import galsim
import mssg_drawLibrary
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- Begin with parameters for the image 

# Parameters for object a (this is a residual from when i was fitting 2 imgs, a & b, with the code this came from)
flux_a = 5e6          # total counts on the image
HLR_a = 3.            # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = 0
y0_a = 0


# Image properties
psf_sigma = 1e-8                          # arcsec
pixel_scale = 0.2                         # arcsec / pixel
imsize = 5*HLR_a/pixel_scale              # Total size of image on a side, in pixels 
dopsfconvln = 'y'
x_len = y_len = imsize

# Galsim function definitions 
galtype = galsim.Gaussian  # Can also be Sersic

# Define a numpy array of params 
param_array = np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])
    
# Random Seed
# Putting in a zero in the arg means it will vary from plot to plot, *not* be exactly same each time 
randomnum = galsim.BaseDeviate(0) 

# ------------------------------------------------------------------------ Create the image 
    
# Create the dataimg -- this calls the galaxy creation routine in this library function.
# Can be chosen to be Gaussian, so needs just flux HLR and standard centroid and ellip params
# Or Sersic, with both expl and bulge comps, see the routine for more details

dataimg = mssg_drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                         imsize,imsize,pixel_scale,galtype,randomnum)     

# Obtain the image bounds and domain information

# An img always has a center and bounds from which lin vecs of ints of
# length the number of pixels in x and y can be made.

# X and Y are meshes (2D mats) made of these ints, in opposite dirs

x_cen,y_cen,x,y,X,Y = mssg_drawLibrary.return_domain(dataimg)
imgarray = dataimg.array # This returns the 2D numpy array assoc with the image - mg
p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])  # These are all pre-defined nums from above, put into a 1D array



# ----------------------------------------------------------------------- Fit
# This declares ann lmfit obj in prep for actually doing the fit
parameters = lmfit.Parameters()

# Add names into the params var, with seed values and enforce min max bounds for the fit values 
parameters.add('flux_a', value=p0[0])   
parameters.add('hlr_a', value=p0[1], min=0.0)
parameters.add('e1_a', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2_a', value=p0[3], min=-1.0, max=1.0)
parameters.add('x0_a',value=p0[4])
parameters.add('y0_a',value=p0[5])

################# Do the actual fit, and assign to a var 'result'

# -- First arg (mssg_drawLibrary.resid_1obj) yields back the MinFunct = function
# that needs to be minimized -- it's the resid between the dataimg
# (which is simulated in our case), and trials using the params
# specified in the params vec

# -- Second arg = params of the MinFunct

# -- Third arg is composed of: (first arg = actual data img, then the
# args that are fed into the MinFunct itself

result = lmfit.minimize(mssg_drawLibrary.resid_1obj,   parameters,   args=(dataimg ,imsize,imsize,pixel_scale, galtype, galtype ))

################# Now print the results to screen
lmfit.report_errors(result.params)

################# Draw a galaxy with the best fit params that came back from the lmfit
best_fit = mssg_drawLibrary.draw_galaxy_1(result.params['flux_a'].value,
                                           result.params['hlr_a'].value,
                                           result.params['e1_a'].value,
                                           result.params['e2_a'].value,
                                           result.params['x0_a'].value,
                                           result.params['y0_a'].value, x_len,y_len, pixel_scale, galtype, dopsfconvln)


############################################ Make some plots of the above

#### Set pltflag
pltflg = 1

if pltflg == 1:
    # Initialize the geometry of the grid
    gs = gridspec.GridSpec(1,3)
    fig = plt.figure(figsize=(10,5))

    # Plot the orig data img
    f1 = fig.add_subplot(gs[0,0])
    plt.title('Binned Image of Galaxy')
    fig1 = f1.imshow(imgarray,interpolation='none',origin='lower')
    plt.colorbar(fig1, shrink=0.5) 

    # Plot the best fit img
    f2 = fig.add_subplot(gs[0,1])
    plt.title('Best Fit')
    fig2 = f2.imshow(best_fit.array,interpolation='none')
    plt.colorbar(fig2, shrink=0.5)

    # Plot the Resid between them
    f3 = fig.add_subplot(gs[0,2])
    plt.title('Resid')
    fig3 = f3.imshow(best_fit.array - imgarray,interpolation='none')
    plt.colorbar(fig3, shrink=0.5)
    plt.show()
    
sys.exit()


