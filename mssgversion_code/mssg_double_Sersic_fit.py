# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 14:55:10 2014

@author: luis
"""

# Based on 2D_lmfit_study_same_params_vary_distance.py by LCA
# Modifications by MSSG
# Start: 8/12/2014

## Trying to fit a single Sersic galaxy

from __future__ import division
import galsim
import mssg_drawLibrary
import matplotlib.pyplot as plt
import lmfit
import numpy as np
import ipdb

# Begin with parameters for the image----------------------------------------

# Parameters for object a
flux_a = 5e4          # total counts on the image
HLR_a = 3.            # arcsec
e1_a = 0.0
e2_a = 0.0


# Image properties
pixel_scale = 0.2     # arcsec / pixel
imsize = 5*HLR_a/pixel_scale            # Total size of image on a side, in pixels - i changed name from size to imsize -mg

# Separation vars - mg 
min_sep = 0.5
max_sep = 4.0
sep_step = 0.5 

sep_array = np.arange(min_sep, max_sep, sep_step) # mg - Changed name from d_coeff
sep_array = sep_array[::-1] # Reverses array

# Galsim function definitions 
# galtype = galsim.Gaussian
galtype = galsim.Sersic # New for Sersic gals - mg

# Arrays for differences
differences = []
errors = []


# Error type to plot
error_types = ['rel_error','abs_error']
error_type = error_types[1]

# Counter to measure where one is at in the loop
count = 0

 # Frequency to make galaxy img plots in this loop (making this one
 # will make a plot at every step, otherwise they will be 'plotfreq'
 # steps apart
plotfreq = 1


# Parameters for object a
flux_a = 5e4          # total counts on the image
HLR_a = 3.            # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = 0
y0_a = 0
        
# Define a numpy array of params - mg
# param_array = np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])    
p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a])
parameters = lmfit.Parameters()

parameters.add('flux_a', value=p0[0])   
parameters.add('hlr_a', value=p0[1], min=0.0)
parameters.add('e1_a', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2_a', value=p0[3], min=-1.0, max=1.0)
parameters.add('x0_a',value=p0[4])
parameters.add('y0_a',value=p0[5])

    #------------------------------------------------------------------------
    # Create the image

    # Random Seed
# Putting in a zero in the arg means it will vary from plot to plot, *not* exactly same each time - mg
dev_1 = galsim.BaseDeviate(0) 

    
# Fill in the images - mg
img1 = mssg_drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        imsize,imsize,pixel_scale,galtype,dev_1)

# Add a Gaussian
galtype = galsim.Gaussian

x0_b = 2*HLR_a

img2 = mssg_drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a, x0_b, y0_a,
                                        imsize,imsize,pixel_scale,galtype,dev_1)

origimg = img1+img2

galtype = galsim.Sersic
    
# Obtain the image bounds and domain information
x_cen,y_cen,x,y,X,Y = mssg_drawLibrary.return_domain(origimg)
H = origimg.array # I believe this returns the 2D numpy array assoc with the image - mg


# Draw the orig img
plt.figure(1)
plt.title('Orig Gal Image')
plt.imshow(H,interpolation='none',origin='lower')
#        plt.plot(x,pos_contour_a,'k',x,neg_contour_a,'k')
#        plt.plot(x,pos_contour_b,'k',x,neg_contour_b,'k')
plt.xlim([1,imsize])
plt.ylim([1,imsize])
plt.show()   


print " ***** Now about to run lmfit "

result = lmfit.minimize(mssg_drawLibrary.resid_1obj,   parameters,   args=(origimg,imsize,imsize,pixel_scale, galtype ))

print " ***** Now about to get best fit params"

best_fit = mssg_drawLibrary.draw_galaxy_1(result.params['flux_a'].value,
                                       result.params['hlr_a'].value,
                                       result.params['e1_a'].value,
                                       result.params['e2_a'].value,
                                       result.params['x0_a'].value,
                                       result.params['y0_a'].value,
                                       imsize,imsize, pixel_scale, galtype)
    
    # Report the parameters to the interpreter screen                        
lmfit.report_errors(result.params)


print " ***** Now about to draw best fit       "

'''  
plt.figure(2)
plt.title('Best Fit')
plt.imshow(best_fit.array,interpolation='none')
plt.show()
'''
  
print " ***** Now about to draw residual  "

resid = best_fit - origimg 

'''
plt.figure(3)
plt.title('Resid')
plt.imshow(resid.array,interpolation='none')
#plt.show()
'''

############################# Draw the 3 plots

fig = plt.figure()

btodsize = 0.5

galpropsStr = '[ '+str(flux_a)+', '+str(HLR_a)+', '+str(btodsize)+', (0, 0) ]'

titleStr = 'Sersic bulge+disk galaxy with parameters: \n [ Photon Flux, HLR (arcsec), bTOd size fraction, (e1, e2) ] = '+ galpropsStr+ '\n Pixel Scale = '+ str(pixel_scale)+' arcsec/pixel; Image size = ' + str(imsize*pixel_scale) + ' arcsec per side '

fig.suptitle(titleStr, fontsize=12)

#Plotting the mixture
ax11 = fig.add_subplot(131)
c1 = ax11.imshow(H, origin='lower')
ax11.set_title('Orig Gal Image')
plt.colorbar(c1, shrink=.5)

#Plotting the fit
ax12 = fig.add_subplot(132)
c2 = ax12.imshow(best_fit.array, origin='lower')
ax12.set_title('Best Fit')
plt.colorbar(c2, shrink=.5)

#Plotting the residual
ax13 = fig.add_subplot(133)
c3 = ax13.imshow((origimg - best_fit).array, origin='lower')
ax13.set_title('Residual')
plt.colorbar(c3, shrink=.5)

plt.show()

# ipdb.set_trace()

