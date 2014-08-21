# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 13:54:38 2014

@author: luis
"""

## This script runs through a simple double fit of two sersic profiles

from __future__ import division
# from pandas import Series, DataFrame
# import pandas as pd
import numpy as np
import mssg_drawLibrary
import mssg_noiseLibrary
import galsim
import lmfit
import matplotlib.pyplot as plt
import scipy.linalg as scla
import matplotlib.gridspec as gridspec
import sys
import ipdb
import time
# import timeit

start= time.time()

# Parameters for object a
flux_a = 5e5          # total counts on the image
hlr_a = 1         # arcsec -- code will segfault without telling you why, if you accidentally set this negative
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 4

# Parameters for object b
flux_b = flux_a       # total counts on the image
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = +1
y0_b = 0
n_b = n_a

# Galsim function definitions
sersic_func = galsim.Sersic
galtype = galsim.Sersic
# galtype = galsim.Gaussian

# Set the RNG
seed_1 = galsim.BaseDeviate(1)
seed_2 = galsim.BaseDeviate(2)
seed_3 = galsim.BaseDeviate(3)

# Image properties
pixel_scale = 0.2     # arcsec / pixel
imsize = 100            # pixels
x_len = y_len = imsize
sky_level = 11         # counts / pixel
add_noise_flag = True


# psf properties
psf_flag = False
beta = 3
fwhm_psf = 0.6


# Separation vars - mg 
min_sep = 0.5
max_sep = 6.0
sep_step = 0.5

sep_array = np.arange(min_sep, max_sep, sep_step) # mg - Changed name from d_coeff

sep_array = sep_array[::-1] # Reverses array

# Arrays for differences
differences = []
errors = []

# Error type to plot
error_types = ['rel_error','abs_error']
error_type = error_types[1]

cur= time.time()

print " About to enter loop -- current-start time = ", cur-start 



for sep in sep_array: 
    loopstarttime = time.time()
    print " At start of sep loop -- current-start time = ", cur-start 
    print " sep = ", sep

    x0_a = -sep*hlr_a  # We'll step gal A to the left in steps of 0.5*step_sep - mg (Shift in arcsec from center)
    y0_a = 0
    x0_b = sep*hlr_b   # We'll step gal B to the right in steps of 0.5*step_sep - mg
    y0_b = 0    

    print " \n\n  Creating Galaxies.. \n\n "
    galimg_a = mssg_noiseLibrary.create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype,sersic_index=n_a, x_len = imsize , y_len = imsize)
    galimg_b = mssg_noiseLibrary.create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype,sersic_index=n_b, x_len = imsize , y_len = imsize)
    
    galimg = galimg_a + galimg_b
    galimg.addNoise(galsim.PoissonNoise(sky_level=0) )

    galimg_array = galimg.array
    

########################################################################################################################################
    
    '''    
# Obtain instantiation
    img, best_fit, result = mssg_noiseLibrary.fit_2_galaxies(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                             flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                             psf_flag,beta,fwhm_psf,
                                                             x_len,y_len,pixel_scale,galtype, galtype, seed_1,seed_2,seed_3,
                                                         add_noise_flag,sky_level, galimg )
    
    
# Get resid
    resid = best_fit - img
    
# Report errors
    lmfit.report_errors(result.params)
    
    print " \n Have fit the 2 gal img \n"
    
    '''

########### Show the galaxies + fit
    fig = plt.figure()

    galpropsStr_a = '[ ('+ str(x0_a) + ',' +str(y0_a) +'), ' + str(n_a) + ', ' + str(flux_a)+', '+str(hlr_a) + ', (0, 0) ]'
    galpropsStr_b = '[ ('+ str(x0_b) + ',' +str(y0_b) +'), '+ str(n_b) + ', ' + str(flux_b)+', '+str(hlr_b) + ', (0, 0) ]'
    
    titleStr = 'Sersic galaxy A with parameters: \n [ (x0,y0) , Sersic Index, Photon Flux, hlr (arcsec), (e1, e2) ] = '+ galpropsStr_a + ' \n And Sersic galaxy B with parameters: \n [  (x0,y0) ,Sersic Index, Photon Flux, hlr (arcsec), (e1, e2) ] = '+ galpropsStr_b + '\n\n Pixel Scale = '+ str(pixel_scale)+' arcsec/pixel; Image size = ' + str(imsize*pixel_scale) + ' arcsec per side '
    
    fig.suptitle(titleStr, fontsize=12)
    
#Plotting the mixture
    ax1 = fig.add_subplot(131)
    c1 = ax1.imshow(galimg_array, origin='lower')
    ax1.set_title('Orig Gals Image')
    plt.colorbar(c1, shrink=.5)

    '''    
    ax2 = fig.add_subplot(132)
    c2 = plt.imshow(best_fit.array, origin='lower')
    ax2.set_title('Best Fit Image')
    plt.colorbar(c2, shrink=.5)
    
    ax3 = fig.add_subplot(133)
    c3 = plt.imshow(resid.array, origin='lower')
    ax3.set_title('Resid Image')
    plt.colorbar(c3, shrink=.5)
   '''
 
    plt.show()
    
    '''
    # Differences on measure and true for object one    
    diff_flux_a = result.params['flux_a'].value - flux_a
    diff_hlr_a = result.params['hlr_a'].value - hlr_a
    diff_e1_a = result.params['e1_a'].value - e1_a
    diff_e2_a = result.params['e2_a'].value - e2_a
    diff_x0_a = result.params['x0_a'].value - x0_a
    diff_y0_a = result.params['y0_a'].value - y0_a
    
    # Differences on measure and true for object two    
    diff_flux_b = result.params['flux_b'].value - flux_b
    diff_hlr_b = result.params['hlr_b'].value - hlr_b
    diff_e1_b = result.params['e1_b'].value - e1_b
    diff_e2_b = result.params['e2_b'].value - e2_b
    diff_x0_b = result.params['x0_b'].value - x0_b
    diff_y0_b = result.params['y0_b'].value - y0_b


    # Append those differences to array
    differences.append([diff_flux_a,diff_hlr_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,
                        diff_flux_b,diff_hlr_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b])

    # Obtain the errors on each parameter
    errors.append(np.sqrt(np.diag(result.covar)))                           
    '''

    cur = time.time()
    print " At end of sep loop -- current-start time = ", cur-loopstarttime 

    print " ********************************* On to next sep.. "

########################## End of loop through separations





# Convert data to numpy arrays for access
errors = np.array(errors)
differences = np.array(differences)
    
# Plot the difference for each estimate on object a and b -------------------

fig = 1
plt.figure(fig)
col = 0
colp = col + 6 # ....? -mg
var_str = '$Flux$'
var = flux_a
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type +' \n $True$ ' + var_str + '$\/= %.2e$' % var)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$Flux_{am}-Flux_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$Flux_{bm}-Flux_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')

fig = 2
plt.figure(fig)
col = 1
colp = col + 6
var_str = '$hlr$'
var = hlr_a
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type +' \n $True$ ' + var_str + '$\/= %.2e$' % var)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$hlr_{am}-hlr_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$hlr_{bm}-hlr_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')
    
fig = 3
plt.figure(fig)
col = 2
colp = col + 6
var_str = '$e_1$'
var = e1_a
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type +' \n $True$ ' + var_str + '$\/= %.2e$' % var)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$e1_{am}-e1_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$e1_{bm}-e1_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')

fig = 4
plt.figure(fig)
col = 3
colp = col + 6
var_str = '$e_2$'
var = e2_a
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type +' \n $True$ ' + var_str + '$\/= %.2e$' % var)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$e2_{am}-e2_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$e2_{bm}-e2_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')
    
fig = 5
plt.figure(fig)
col = 4
colp = col + 6
var_str = '$x0$'
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$x0_{am}-x0_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$x0_{bm}-x0_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')
    
fig = 6
plt.figure(fig)
col = 5
colp = col + 6
var_str = '$y0$'
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$y0_{am}-y0_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$y0_{bm}-y0_{bt}$')
plt.errorbar(sep_array,differences[:,colp],yerr=errors[:,colp])    
plt.legend()
plt.xlabel('$Distance\/Coefficent\/Between\/Galaxies\/$')
plt.ylabel('$Difference$')
    
plt.show()
    

sys.exit()
################################################################################################################

im = img

# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(1,3)
# Set the figure size
fig = plt.figure(figsize=(15,10))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+10)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,0])
plt.title('$Image$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im.array,interpolation='none',origin='lower',vmax=im.array.max(),vmin=im.array.min())
plt.colorbar(a,shrink=0.5)
# Add the second subplot
ax2 = fig.add_subplot(gs[0,1])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Best Fit$',fontsize=fonts+8)
b = ax2.imshow(best_fit.array,origin='lower',vmax=best_fit.array.max(),vmin=best_fit.array.min())
plt.colorbar(b,shrink=0.5)
# Add the third subplot
ax3 = fig.add_subplot(gs[0,2])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Residual$',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow((best_fit - im).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=0.5)

ipdb.set_trace()                              
