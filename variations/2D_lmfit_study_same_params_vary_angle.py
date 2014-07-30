# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 10:00:16 2014

@author: luis
"""

## Loop through different angles for a secondary galaxy
## with the exact profile as the same as the primary galaxy.

from __future__ import division
import galsim
import drawLibrary
import matplotlib.pyplot as plt
import lmfit
import numpy as np

# Angle to create the secondary galaxy at with respect to the center
# of the image.
thetas = np.linspace(0,360-22.5,16)

# Begin with parameters for the image

# Parameters for object a
flux_a = 5e6       # total counts on the image
HLR_a = 3.          # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = 0
y0_a = 0

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 1/5     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
size = 500             # pixel

# Galsim function definitions
func_gauss = galsim.Gaussian

# Count 
count = 0

# Arrays for differences
differences = []
errors = []

# Arrays for different angles
ninety_deg_axes = []; error_ninety = []
forty_five_deg_axes = []; error_forty_five = []
twenty_two_five_deg_axes = []; error_twenty_two_five = []

# Distance factor
d = 1

# Frequency of plots for loop
freq = 18

# Error type to plot
error_types = ['rel_error','abs_error']
error_type = error_types[1]

# Counter to measure where one is at in the loop
count = 0

for theta in thetas: 
    
    org_theta = theta
    theta += 0.0

    # Parameters for object b
    flux_b = flux_a       # total counts on the image
    HLR_b = HLR_a          # arcsec
    e1_b = 0.0
    e2_b = 0.0
    x0_b = d*(HLR_a + HLR_b)*np.cos(theta*np.pi/180)       # Shift of arcsec from center
    y0_b = d*(HLR_a + HLR_b)*np.sin(theta*np.pi/180)        # Shift of arcsec from center
    
    param_array = np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                            flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b])    
    
    #------------------------------------------------------------------------
    # Create the image
    
    # Random Seed
    dev_1 = galsim.BaseDeviate(0)
    dev_2 = galsim.BaseDeviate(0)
    
    im_1 = drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        size,size,pixel_scale,func_gauss,dev_1)
    im_2 = drawLibrary.drawShoot_galaxy(flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                        size,size,pixel_scale,func_gauss,dev_2)
    im = im_1 + im_2                                
    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(im)
    H = im.array

    # Create the contours of the two galaxies    
    pos_contour_a = drawLibrary.contour_positive_circle(x,x_cen,y_cen,HLR_a/pixel_scale)
    neg_contour_a = drawLibrary.contour_negative_circle(x,x_cen,y_cen,HLR_a/pixel_scale)
 
    pos_contour_b = drawLibrary.contour_positive_circle(x,x_cen+x0_b/pixel_scale,y_cen+y0_b/pixel_scale,HLR_b/pixel_scale)
    neg_contour_b = drawLibrary.contour_negative_circle(x,x_cen+x0_b/pixel_scale,y_cen+y0_b/pixel_scale,HLR_b/pixel_scale)           
    
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
          flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
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
    result = lmfit.minimize(drawLibrary.resid_2, parameters, args=(im,size,size,pixel_scale,func_gauss,func_gauss))
    best_fit = drawLibrary.draw_galaxy_2(result.params['flux_1'].value,
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
                                       size,size,pixel_scale,func_gauss,func_gauss)
    
    # Report the parameters to the interpreter screen                        
    lmfit.report_errors(result.params)
    
    # Print the 
    print 'Angle for plot = '+ str(org_theta)
    print 'Angle in calculation = ' + str(theta)
    
    # Update the count
    count += 1
    
    # Intermediate plots to display galaxies and their best fit
    if np.mod(count,freq) == 0:
        # Plot the binned data as looping occurs        
        plt.figure(1)
        plt.title('Binned Image of Galaxies')
        plt.imshow(H,interpolation='none',origin='lower')
        plt.plot(x,pos_contour_a,'k',x,neg_contour_a,'k')
        plt.plot(x,pos_contour_b,'k',x,neg_contour_b,'k')
        plt.xlim([1,size])
        plt.ylim([1,size])
        plt.show()   
        # Plot the best fit as looping occur
        plt.figure(2)
        plt.title('Best Fit')
        plt.imshow(best_fit.array,interpolation='none',origin='lower')
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

    # Use algorithm to sort entries by pure multiple of pi/8,pi/4,pi/2 ------
    
    # If the count is odd, insert into 0 or pi/4 category
    if np.mod(count,2) == 1:
        # If theta is a multiple of 90.1 then insert into 0 axes
        if np.floor(np.mod(org_theta,90)) == 0 or theta == 0:
            ninety_deg_axes.append([diff_flux_a,diff_HLR_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,
                            diff_flux_b,diff_HLR_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b])
            error_ninety.append(np.sqrt(np.diag(result.covar)))               
                   
        # If theta is a multiple of 45, insert into other        
        elif np.mod(org_theta,45) == 0:
            forty_five_deg_axes.append([diff_flux_a,diff_HLR_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,
                            diff_flux_b,diff_HLR_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b])                    
            error_forty_five.append(np.sqrt(np.diag(result.covar)))
    # If the count is even, then theta is a multiple of 22.5
    elif np.mod(count,2) == 0:
        twenty_two_five_deg_axes.append([diff_flux_a,diff_HLR_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,
                            diff_flux_b,diff_HLR_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b])
        error_twenty_two_five.append(np.sqrt(np.diag(result.covar)))                    

# Convert data to numpy arrays for access
errors = np.array(errors)
differences = np.array(differences)
    
# Plot the difference for each estimate on object a and b -------------------

plt.figure(1)
plt.suptitle('Difference of $Flux$ for Distance Factor of $%d$. Error Type: %s\n $True\/Flux = %s$' % (d,error_type,flux_a))
plt.scatter(thetas,differences[:,0],c='b',marker='x',label='$Flux_{am}-Flux_{at}$')
plt.errorbar(thetas,differences[:,0],yerr=errors[:,0])  
plt.scatter(thetas,differences[:,6],c='g',marker='o',label='$Flux_{bm}-Flux_{bt}$')
plt.errorbar(thetas,differences[:,6],yerr=errors[:,6])    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')
    

plt.figure(2)
plt.suptitle('Difference of $HLR$ for Distance Factor of $%d$. Error Type: %s\n $True\/HLR = %s$' % (d,error_type,HLR_a))
plt.scatter(thetas,differences[:,1],c='b',marker='x',label='$HLR_{am}-HRL_{at}$')
plt.errorbar(thetas,differences[:,1],yerr=errors[:,1],ecolor='blue')  
plt.scatter(thetas,differences[:,7],c='g',marker='o',label='$HLR_{bm}-HRL_{bt}$')
plt.errorbar(thetas,differences[:,7],yerr=errors[:,7],ecolor='green')    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')

plt.figure(3)
plt.suptitle('Difference of $e1$ for Distance Factor of $%d$. Error Type: %s\n $True\/e_1 = %s$' % (d,error_type,e1_a))
plt.scatter(thetas,differences[:,2],c='b',marker='x',label='$e1_{am}-e1_{at}$')
plt.errorbar(thetas,differences[:,2],yerr=errors[:,2])  
plt.scatter(thetas,differences[:,8],c='g',marker='o',label='$e1_{bm}-e1_{bt}$')
plt.errorbar(thetas,differences[:,8],yerr=errors[:,8])    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')

plt.figure(4)
plt.suptitle('Difference of $e2$ for Distance Factor of $%d$. Error Type: %s\n $True\/e_2 = %s$' % (d,error_type,e2_a))
plt.scatter(thetas,differences[:,3],c='b',marker='x',label='$e2_{am}-e2_{at}$')
plt.errorbar(thetas,differences[:,3],yerr=errors[:,3])  
plt.scatter(thetas,differences[:,9],c='g',marker='o',label='$e2_{bm}-e2_{bt}$')
plt.errorbar(thetas,differences[:,9],yerr=errors[:,9])    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')
        
plt.figure(5)
plt.suptitle('Difference of $x0$ for Distance Factor of $%d$. Error Type: %s\n' % (d,error_type))
plt.scatter(thetas,differences[:,4],c='b',marker='x',label='$x0_{am}-x0_{at}$')
plt.errorbar(thetas,differences[:,4],yerr=errors[:,4])  
plt.scatter(thetas,differences[:,10],c='g',marker='o',label='$x0_{bm}-x0_{bt}$')
plt.errorbar(thetas,differences[:,10],yerr=errors[:,10])    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')

plt.figure(6)
plt.suptitle('Difference of $y0$ for Distance Factor of $%d$. Error Type: %s' % (d,error_type))
plt.scatter(thetas,differences[:,5],c='b',marker='x',label='$y0_{am}-y0_{at}$')
plt.errorbar(thetas,differences[:,5],yerr=errors[:,5])  
plt.scatter(thetas,differences[:,11],c='g',marker='o',label='$y0_{bm}-y0_{bt}$')
plt.errorbar(thetas,differences[:,11],yerr=errors[:,11])    
plt.legend()
plt.xlabel('$Angle\/of\/Secondary\/Galaxy\/(deg)$')
plt.ylabel('$Difference$')
plt.show()

# Plot the new data with both the average and all the points --------------

# Convert the data arrays to numpy arrays for indexing
twenty_two_five_deg_axes = np.array(twenty_two_five_deg_axes)
forty_five_deg_axes = np.array(forty_five_deg_axes)
ninety_deg_axes = np.array(ninety_deg_axes)

# Plot the Flux difference ----------------
plt.figure(7)
col = 0 # Column number to choose data from
var_str = '$Flux$'
offset = 0.05
length_col = 4
var = flux_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%.1e$'%(var))
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)

# Plot the HLR ----------------------------

plt.figure(8)
col = 1 # Column number to choose data from
var_str = '$HLR_a$'
offset = 0.05
length_col = 4
var = HLR_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%d$'%(var))
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)


# Plot e1 ---------------------------------

plt.figure(9)
col = 2 # Column number to choose data from
var_str = '$e1_a$'
offset = 0.05
length_col = 4
var = e1_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%d$'%(var))
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)

# Plot e2 ---------------------------------

plt.figure(10)
col = 3 # Column number to choose data from
var_str = '$e2_a$'
offset = 0.05
length_col = 4
var = e2_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%d$'%(var))
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)

# Plot x0 ---------------------------------

plt.figure(11)
col = 4 # Column number to choose data from
var_str = '$x0_a$'
offset = 0.05
length_col = 4
var = x0_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%d$'%(var))
# 22.5 Degrees
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)

# Plot y0 ---------------------------------

plt.figure(12)
col = 5 # Column number to choose data from
var_str = '$y0_a$'
offset = 0.05
length_col = 4
var = y0_a
plt.suptitle('Difference of ' + var_str +' vs Angle Bin\n True ' + var_str + ' $=\/%d$'%(var))
# 22.5 Degrees
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][:4],marker='s',c='b')
plt.scatter(np.ones(length_col),twenty_two_five_deg_axes[:,col][4:8],marker='o',c='b')
#plt.errorbar(np.ones(2*length_col),twenty_two_five_deg_axes[:,col],yerr=twenty_two_five_deg_axes[:,col],linestyle='none')
# 45 Degrees
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][:2],marker='s',c='g')
plt.scatter(2*np.ones(length_col/2),forty_five_deg_axes[:,col][2:4],marker='o',c='g')
#plt.errorbar(2*np.ones(length_col),forty_five_deg_axes[:,col],yerr=forty_five_deg_axes[:,col],linestyle='none')
# 90 Degrees
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][:2],marker='s',c='r')
plt.scatter(3*np.ones(length_col/2),ninety_deg_axes[:,col][2:4],marker='o',c='r')
#plt.errorbar(3*np.ones(length_col),ninety_deg_axes[:,col],yerr=ninety_deg_axes[:,col],linestyle='none')

# 22.5 Degrees
p1 = plt.scatter(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(np.ones(1),np.mean(twenty_two_five_deg_axes[:,col]),yerr=np.std(twenty_two_five_deg_axes[:,col]),linestyle='none')
# 45 Degrees
p2 = plt.scatter(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(2*np.ones(1),np.mean(forty_five_deg_axes[:,col]),yerr=np.std(forty_five_deg_axes[:,col]),linestyle='none')
# 90 Degrees
p3 = plt.scatter(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),marker='x',c='k',linewidths=6)
plt.errorbar(3*np.ones(1),np.mean(ninety_deg_axes[:,col]),yerr=np.std(ninety_deg_axes[:,col]),linestyle='none')
plt.text(1+offset,np.mean(twenty_two_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(2+offset,np.mean(forty_five_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
plt.text(3+offset,np.mean(ninety_deg_axes[:,col]),'Mean Value',fontsize=10,zorder=10)
#plt.legend([p1,p2,p3],['$Mean\/Value$','$Mean\/Value$','$Mean\/Value$'])

plt.ylabel('$Difference$')
plt.xticks([1.0,2.0,3.0],[r'$\frac{\pi}{8}$',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$'],fontsize=15)
plt.xlabel('$Angle\/Bin$',fontsize=11)

plt.show()