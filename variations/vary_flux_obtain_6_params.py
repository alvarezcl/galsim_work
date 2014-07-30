# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 16:25:02 2014

@author: luis
"""

from __future__ import division
import galsim
import drawLibrary
import matplotlib.pyplot as plt
import lmfit
import numpy as np
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

# Function definition to return the original data array, best-fit array,
# residual, and correlation matrix with differences and error on e1 and e2.
def run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                               flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                               size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2,dev_1,dev_2):

    im_1 = drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        size,size,pixel_scale,func_gauss,dev_1)
    im_2 = drawLibrary.drawShoot_galaxy(flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                        size,size,pixel_scale,func_gauss,dev_2)
    
    im = im_1 + im_2
    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(im)
    
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
    result = lmfit.minimize(drawLibrary.resid_2, parameters, args=(im,size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2))
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
                                       size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2)
    
#    if result.covar is None:
#        print "Undefined Covariance Matrix\n Rerunning"
#        return run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
#                                          flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
#                                          size_1,size_2,pixel_scale,
#                                          func_gauss_1,func_gauss_2,galsim.BaseDeviate(0),galsim.BaseDeviate(0))
#    else:
#        print "result.covar is defined."    
    
#    error_diag = np.sqrt(np.diag(result.covar))
#    error_mat = np.outer(error_diag,error_diag)
#    correlation_mat = result.covar/error_mat
    diff_flux_a = result.params['flux_1'].value - flux_a
    diff_HLR_a = result.params['hlr_1'].value - HLR_a
    diff_e1_a = result.params['e1_1'].value - e1_a
    diff_e2_a = result.params['e2_1'].value - e2_a
    diff_x0_a = result.params['x_center1'].value - x0_a
    diff_y0_a = result.params['y_center1'].value - y0_a 

    diff_flux_b = result.params['flux_2'].value - flux_b
    diff_HLR_b = result.params['hlr_2'].value - HLR_b                              
    diff_e1_b = result.params['e1_2'].value - e1_b
    diff_e2_b = result.params['e2_2'].value - e2_b
    diff_x0_b = result.params['x_center2'].value - x0_b
    diff_y0_b = result.params['y_center2'].value - y0_b 

    
    
    diff = (diff_flux_a,diff_HLR_a,diff_e1_a,diff_e2_a,diff_x0_a,diff_y0_a,diff_flux_b,diff_HLR_b,diff_e1_b,diff_e2_b,diff_x0_b,diff_y0_b)
    #error_e = [error_diag[2],error_diag[3],error_diag[8],error_diag[9]]
    
    #error_e = np.array(error_e)
                                                                      
    return im.array,best_fit.array,(best_fit-im).array,diff,result.covar


# Begin with parameters for the image----------------------------------------

# Parameters for object a
flux_a = 5e6          # total counts on the image
HLR_a = 1.            # arcsec
e1_a = 0.0
e2_a = 0.0

# Parameters for object b
flux_b = flux_a       # total counts on the image
HLR_b = HLR_a         # arcsec
e1_b = 0.0
e2_b = 0.0

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 1/5     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
size = 100            # pixel

# Set the beginning locations
x0_a = -HLR_a
y0_a = 0

x0_b = HLR_b
y0_b = 0.0

distance = x0_b - x0_a # arcsec

# Galsim function definitions
func_gauss = galsim.Gaussian

# Set the RNG
dev_1 = galsim.UniformDeviate(1)
dev_2 = galsim.UniformDeviate(2)

# Num of frames from animation
frames = 97

# Separation array
flux_r = []
resid = []

# Run through loop to get residuals
for i in xrange(0,frames):
    
    flux_ratio = flux_b/flux_a 
    flux_r.append(flux_ratio)
    print flux_ratio
    im, best_fit, residual,diff,result_covar = run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                                             flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                                             size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2)                                                                                    
    while result_covar is None:
        print "Undefined Covariance Matrix\n Rerunning"
        im, best_fit, residual,diff,result_covar = run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                                             flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                                             size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2)
    else:
        print "result.covar is defined."    
    
    error_diag = np.sqrt(np.diag(result_covar))
    resid.append([diff[0],diff[1],diff[2],diff[3],diff[4],diff[5],diff[6],diff[7],diff[8],diff[9],diff[10],diff[11]])
        
    flux_b += -50000

flux_r = np.array(flux_r)
resid = np.array(resid)
    
# Provide the fontsize
fonts = 20
# Initialize the geometry of the grid
gs = gridspec.GridSpec(9,1)
# Figure
fig = plt.figure(figsize=(20,11))
# Column to choose from
col = 3
col = col - 1
colp = col + 6
# Figure
ax1 = fig.add_subplot(gs[:4,0])    
plt.title('Residual of $e_1$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlim([0,1.2])
#plt.xlabel('$Separation\/(Arcsec)$ ',fontsize=fonts)
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$e1_{A}^{(m)}-e1_{A}^{(t)}$','$e1_{B}^{(m)}-e1_{B}^{(t)}$'],loc=1,prop={'size':fonts/1.5})

# Column to choose from
col = 4
col = col - 1
colp = col + 6
# Figure    
ax2 = fig.add_subplot(gs[5:,0])
plt.title('Residual of $e_2$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlim([0,1.2])
plt.xlabel('$Flux\/Ratio\/$ ',fontsize=fonts)
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$e2_{A}^{(m)}-e2_{A}^{(t)}$','$e2_{B}^{(m)}-e2_{B}^{(t)}$'],loc=1,prop={'size':fonts/1.5})
plt.show()


# Figure
fig1 = plt.figure(figsize=(20,11))
# Column to choose from
col = 1
col = col - 1
colp = col + 6
# Figure    
ax2 = fig1.add_subplot(gs[:4,0])
plt.title('Residual of $Flux$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlim([0,1.2])
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$Flux_{A}^{(m)}-Flux_{A}^{(t)}$','$Flux_{B}^{(m)}-Flux_{B}^{(t)}$'],loc=1,prop={'size':fonts/1.5})

# Column to choose from
col = 2
col = col - 1
colp = col + 6
# Figure
ax1 = fig1.add_subplot(gs[5:,0])    
plt.title('Residual of $HLR$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlim([0,1.2])
plt.xlabel('$Flux\/Ratio\/$ ',fontsize=fonts)
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$HLR_{A}^{(m)}-HLR_{A}^{(t)}$','$HLR_{B}^{(m)}-HLR_{B}^{(t)}$'],loc=4,prop={'size':fonts/1.5})
plt.show()

# Figure
fig2 = plt.figure(figsize=(20,11))
# Column to choose from
col = 5
col = col - 1
colp = col + 6
# Figure
ax1 = fig2.add_subplot(gs[:4,0])    
plt.title('Residual of $x0$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlim([0,1.2])
#plt.xlabel('$Separation\/(Arcsec)$ ',fontsize=fonts)
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$x0_{A}^{(m)}-x0_{A}^{(t)}$','$x0_{B}^{(m)}-y0_{B}^{(t)}$'],loc=1,prop={'size':fonts/1.5})

# Column to choose from
col = 6
col = col - 1
colp = col + 6
# Figure    
ax2 = fig2.add_subplot(gs[5:,0])
plt.xlim([0,1.2])
plt.title('Residual of $y0$ vs $Flux\/Ratio$ Of Object $B$',fontsize=fonts)
plt.xlabel('$Flux\/Ratio\/$ ',fontsize=fonts)
plt.ylabel('$Residual$',fontsize=fonts)
p1 = plt.plot(flux_r,resid[:,col],c='b')
p2 = plt.plot(flux_r,resid[:,colp],c='g')
p1 = plt.scatter(flux_r,resid[:,col],c='b',marker='o')
p2 = plt.scatter(flux_r,resid[:,colp],c='g',marker='x')
plt.legend([p1,p2],['$y0_{A}^{(m)}-y0_{A}^{(t)}$','$y0_{B}^{(m)}-y0_{B}^{(t)}$'],loc=4,prop={'size':fonts/1.5})
plt.show()