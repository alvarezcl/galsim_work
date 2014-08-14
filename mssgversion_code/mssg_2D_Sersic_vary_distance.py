# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 14:55:10 2014

@author: luis
"""

# Based on 2D_lmfit_study_same_params_vary_distance.py by LCA
# Modifications by MSSG
# Start: 8/2/2014

## Loop through different distances for a secondary galaxy
## with profile same as for the primary galaxy.

from __future__ import division
import galsim
import mssg_drawLibrary
import matplotlib.pyplot as plt
import lmfit
import numpy as np
import ipdb

# Begin with parameters for the image----------------------------------------

# Parameters for object a
flux_a = 5e6          # total counts on the image
HLR_a = 3.            # arcsec
e1_a = 0.0
e2_a = 0.0

# Parameters for object b
flux_b = flux_a       # total counts on the image
HLR_b = HLR_a         # arcsec
e1_b = 0.0
e2_b = 0.0

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 0.2     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
imsize = 5*HLR_a/pixel_scale            # Total size of image on a side, in pixels - i changed name from size to imsize -mg

# Separation vars - mg 
min_sep = 0.5
max_sep = 6.0
sep_step = 0.5 # This makes it 16 steps

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

###############  Step through the separations in the sep_array -mg
for sep in sep_array: 
    
    print " sep = ", sep

    x0_a = -sep*HLR_a*0.5  # We'll step gal A to the left in steps of 0.5*step_sep - mg (Shift in arcsec from center)
    y0_a = 0
    x0_b = sep*HLR_b*0.5   # We'll step gal B to the right in steps of 0.5*step_sep - mg
    y0_b = 0    

    
         
# Define a numpy array of params - mg
    param_array = np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                            flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b])    
    
    #------------------------------------------------------------------------
    # Create the image
    
    # Random Seed
# Putting in a zero in the arg means it will vary from plot to plot, *not* exactly same each time - mg
    dev_1 = galsim.BaseDeviate(0) 
    dev_2 = galsim.BaseDeviate(0)    
    
# Fill in the images - mg
    im_1 = mssg_drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        imsize,imsize,pixel_scale,galtype,dev_1)
    im_2 = mssg_drawLibrary.drawShoot_galaxy(flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                        imsize,imsize,pixel_scale,galtype,dev_2)
    
    im = im_1 + im_2    
    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = mssg_drawLibrary.return_domain(im)
    H = im.array # I believe this returns the 2D numpy array assoc with the image - mg



    # Create the contours of the two galaxies    
#    pos_contour_a = mssg_drawLibrary.contour_positive_circle(x,x_cen+x0_a/pixel_scale,y_cen+y0_a/pixel_scale,HLR_a/pixel_scale)
#    neg_contour_a = mssg_drawLibrary.contour_negative_circle(x,x_cen+x0_a/pixel_scale,y_cen+y0_a/pixel_scale,HLR_a/pixel_scale)
 
#    pos_contour_b = mssg_drawLibrary.contour_positive_circle(x,x_cen+x0_b/pixel_scale,y_cen+y0_b/pixel_scale,HLR_b/pixel_scale)
#    neg_contour_b = mssg_drawLibrary.contour_negative_circle(x,x_cen+x0_b/pixel_scale,y_cen+y0_b/pixel_scale,HLR_b/pixel_scale)           

    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
          flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b]) # These are all pre-defined nums from above - mg
    parameters = lmfit.Parameters()
# Shouldn't these really be all called _a and then _b vs. 1 and 2..?  Maybe they're just name labels though..?  -mg
    parameters.add('flux_1', value=p0[0])   
    parameters.add('hlr_1', value=p0[1], min=0.0)
    parameters.add('e1_1', value=p0[2], min=-1.0, max=1.0)
    parameters.add('e2_1', value=p0[3], min=-1.0, max=1.0)
    parameters.add('x_center1',value=p0[4])
    parameters.add('y_center1',value=p0[5])
    
    parameters.add('flux_2', value=p0[6])
    parameters.add('hlr_2', value=p0[7], min=0.0)
    parameters.add('e1_2', value=p0[8], min=-1.0, max=1.0)
    parameters.add('e2_2', value=p0[9], min=-1.0, max=1.0)
    parameters.add('x_center2',value=p0[10])
    parameters.add('y_center2',value=p0[11])
    
    
    # Extract params that minimize the difference of the data from the model.

# This is how lmfit is called: 
# -- the first arg is the thing to be minimized
# -- second is the list of params, 
# -- third is the list of args we pass to it, which are the image and its properties, including the function we're fitting to for each obj

    print " About to run lmfit "

### Set the type of galaxy to fit
#    galtype = galsim.Sersic 

    result = lmfit.minimize(mssg_drawLibrary.resid_2obj,   parameters,   args=(im,imsize,imsize,pixel_scale, galtype, galtype ))


    best_fit = mssg_drawLibrary.draw_2galaxies(result.params['flux_1'].value,
                                       result.params['hlr_1'].value,
                                       result.params['e1_1'].value,
                                       result.params['e2_1'].value,
                                       result.params['x_center1'].value,
                                       result.params['y_center1'].value,
                                       result.params['flux_2'].value,
                                       result.params['hlr_2'].value,
                                       result.params['e1_2'].value,
                                       result.params['e2_2'].value,
                                       result.params['x_center2'].value,
                                       result.params['y_center2'].value,
                                       imsize,imsize,pixel_scale, galtype, galtype)
    
    # Report the parameters to the interpreter screen                        
    lmfit.report_errors(result.params)
    
    # Update the loop count, so we know whether to make plots in the next if-then -mg
    count += 1
    
    # Intermediate plots to display galaxies and their best fit
    if np.mod(count,plotfreq) == 0:
        # Plot the binned data as looping occurs        
        plt.figure(1)
        plt.title('Binned Image of Galaxies')
        plt.imshow(H,interpolation='none',origin='lower')
#        plt.plot(x,pos_contour_a,'k',x,neg_contour_a,'k')
#        plt.plot(x,pos_contour_b,'k',x,neg_contour_b,'k')
        plt.xlim([1,imsize])
        plt.ylim([1,imsize])
        plt.show()   
        # Plot the best fit as looping occur
        plt.figure(2)
        plt.title('Best Fit')
        plt.imshow(best_fit.array,interpolation='none')
        plt.show()
        error_diag = np.sqrt(np.diag(result.covar))
        error_mat = np.outer(error_diag,error_diag)
        correlation_mat = result.covar/error_mat
        plt.title('Correlation Coefficient Matrix')
        plt.imshow(correlation_mat,interpolation='none',origin='lower',vmin=-1,vmax=1)
        plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
                   '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'])
        plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
                   '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'])
        plt.colorbar()
        plt.show()
    
    # Differences on measure and true for object one    
    diff_flux_a = result.params['flux_1'].value - flux_a
    diff_HLR_a = result.params['hlr_1'].value - HLR_a
    diff_e1_a = result.params['e1_1'].value - e1_a
    diff_e2_a = result.params['e2_1'].value - e2_a
    diff_x0_a = result.params['x_center1'].value - x0_a
    diff_y0_a = result.params['y_center1'].value - y0_a
    
    # Differences on measure and true for object two    
    diff_flux_b = result.params['flux_2'].value - flux_b
    diff_HLR_b = result.params['hlr_2'].value - HLR_b
    diff_e1_b = result.params['e1_2'].value - e1_b
    diff_e2_b = result.params['e2_2'].value - e2_b
    diff_x0_b = result.params['x_center2'].value - x0_b
    diff_y0_b = result.params['y_center2'].value - y0_b


    # Append those differences to array
    differences.append([diff_flux_a,diff_HLR_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,
                        diff_flux_b,diff_HLR_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b])

    # Obtain the errors on each parameter
    errors.append(np.sqrt(np.diag(result.covar)))                           

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
var_str = '$HLR$'
var = HLR_a
plt.suptitle('Difference of ' + var_str + ' for Changing Distance. Error Type: '+ error_type +' \n $True$ ' + var_str + '$\/= %.2e$' % var)
plt.scatter(sep_array,differences[:,col],c='b',marker='x',label='$HLR_{am}-HLR_{at}$')
plt.errorbar(sep_array,differences[:,col],yerr=errors[:,col])  
plt.scatter(sep_array,differences[:,colp],c='g',marker='o',label='$HLR_{bm}-HLR_{bt}$')
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
    
    
    
    
    
