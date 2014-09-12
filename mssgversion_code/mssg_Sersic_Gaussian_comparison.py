# 8/24/2014
# MSSG

# Code to compare a Sersic vs. a Gaussian profile

from __future__ import division
# from pandas import Series, DataFrame
# import pandas as pd
import numpy as np
import drawLibrary
#import mssg_noiseLibrary as noiseLibrary
import noiseLibrary 
import galsim
import lmfit
import matplotlib.pyplot as plt
import scipy.linalg as scla
import matplotlib.gridspec as gridspec
import sys
import ipdb
import time
import gauss
import plotutils
from scipy.optimize import curve_fit

start= time.time()

# Parameters for object a
flux_a = 5e5          # total counts on the image
hlr_a = 1         # arcsec -- code will segfault without telling you why, if you accidentally set this negative
e1_a = 0.0
e2_a = 0.0
x0_a = 0
y0_a = 0
n_a = 0.5   # 0.5 = Gaussian

# Parameters for object b
flux_b = flux_a       # total counts on the image
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = -x0_a
y0_b = 0
n_b = 1     # 1 = expl disk, 4 = bulge

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
imsize = 50            # pixels
x_len = y_len = imsize
sky_level = 0         # counts / pixel
add_noise_flag = False

cur= time.time()

############################ Plot analytic functions
'''
fig = plt.figure()


x = np.arange(0,10.0, 0.2)

if (n_b == 1 ):
    plotlabel = 'Expl'
    sers = x*np.exp(-x)  # Correct nzn for expl
if (n_b == 4 ):
    plotlabel = 'DeVauc'
    sers = x*np.exp(-(x**0.25))  

gaus = (2/np.sqrt(np.pi))*x*np.exp(-(x**2)) # --> integ to x=0.48  gives 0.504 as area under this curve
plt.plot(x, gaus, label='Gaussian' )
plt.plot(x, sers, label= plotlabel )
plt.legend(loc='upper right')

plt.show()

sys.exit()
'''

############################ Make gals
print " \n\n  Creating Galaxies.. \n\n "
galimg_a = noiseLibrary.create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype,sersic_index=n_a, x_len = imsize , y_len = imsize)
galimg_b = noiseLibrary.create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype,sersic_index=n_b, x_len = imsize , y_len = imsize)
    
difimg = galimg_a - galimg_b

# Make arrays
galimg_a_array = galimg_a.array
galimg_b_array = galimg_b.array


############################ Pick slices
fig = plt.figure() # Have to remake a new fig object

x_cen,y_cen,xx,yy,X,Y = drawLibrary.return_domain(galimg_a)

slicerow  = 25  # 25 = center of img if size is 50 pixels on a side
print " y_cen = ", y_cen 
slice_a = galimg_a.array[y_cen] 
slice_b = galimg_b.array[y_cen] 

xr = xx-x_cen

rad_a = galimg_a.array[slicerow] *xr
rad_b = galimg_b.array[slicerow] *xr


if (n_b == 1 ):
    plotlabel = 'Expl Slice'
if (n_b == 4 ):
    plotlabel = 'DeVauc Slice'


#ipdb.set_trace()

# Plot slices
sp1 = fig.add_subplot(121)
'''
sp1.bar(range(0,imsize), slice_a,alpha=0.5, label='Gaussian Slice' ,color='b')
sp1.bar(range(0,imsize), slice_b,alpha=0.5, label=plotlabel, color='g')
'''
plt.plot(xr, slice_a, alpha=0.5, label='Gaussian Slice' ,color='b')
plt.plot(xr, slice_a*xr, alpha=0.5, label='Gaussian Slice*xr' ,color='g')

plt.plot(xr, slice_b, alpha=0.5, label='Expl Slice' ,color='r')
plt.plot(xr, slice_b*xr, alpha=0.5, label='Expl Slice*xr' ,color='k')

sp1.legend(loc='upper right')

############################ Make projections

proj_a = galimg_a_array.sum(axis=0)
proj_b = galimg_b_array.sum(axis=0)

# Plot projs
sp2 = fig.add_subplot(122)
#sp2.bar(range(0,imsize), proj_a,alpha=0.5, label='Gaussian Proj ' ,color='b')
#sp2.bar(range(0,imsize), proj_b,alpha=0.5, label=plotlabel, color='g')

sp2.bar(range(0,imsize), rad_a,alpha=0.5, label='Gaussian Proj ' ,color='b')
sp2.bar(range(0,imsize), rad_b,alpha=0.5, label='Expl Proj', color='g')

sp2.legend(loc='upper right')

plt.show()

# sys.exit()
#ipdb.set_trace()

########### Show the galaxies + fit
fig = plt.figure()

galpropsStr_a = '[ ('+ str(x0_a) + ',' +str(y0_a) +'), ' + str(n_a) + ', ' + str(flux_a)+', '+str(hlr_a) + ', (0, 0) ]'

galpropsStr_b = '[ ('+ str(x0_b) + ',' +str(y0_b) +'), '+ str(n_b) + ', ' + str(flux_b)+', '+str(hlr_b) + ', (0, 0) ]'
    
# Make title
titleStr = 'Sersic galaxy A with parameters: \n [ (x0,y0) , Sersic Index, Photon Flux, hlr (arcsec), (e1, e2) ] = '+ galpropsStr_a + ' \n\n And Sersic galaxy B with parameters: \n [  (x0,y0) ,Sersic Index, Photon Flux, hlr (arcsec), (e1, e2) ] = '+ galpropsStr_b + '\n\n Pixel Scale = '+ str(pixel_scale)+' arcsec/pixel; Image size = ' + str(imsize*pixel_scale) + ' arcsec per side '
    
fig.suptitle(titleStr, fontsize=12)
    
# Needed to set the color bars for the 2 dif gal profile plots
sersprofilemax = np.max(galimg_b.array)  

#Plot Gaussian
sp1 = fig.add_subplot(131)
c1 = sp1.imshow(galimg_a.array, origin='lower', vmin=0,vmax=sersprofilemax)
sp1.set_title('Galaxy_A Image')
plt.colorbar(c1, shrink=.5)

#Plot Sersic
sp2 = fig.add_subplot(132)
c2 = plt.imshow(galimg_b.array, origin='lower',vmin=0,vmax=sersprofilemax)
sp2.set_title('Galaxy_B Image')
plt.colorbar(c2, shrink=.5)

#Plot dif (A - B)
sp3 = fig.add_subplot(133)
c3 = plt.imshow(difimg.array, origin='lower')
sp3.set_title('A-B Image')
plt.colorbar(c3, shrink=.5)

plt.show()

'''
#### 1D case -- Beginning from an analytic expression, plot histogram

# Input vars
mean = 0; var = 1; sigma = np.sqrt(var); N = 10000; nbins = N/100


# Note we choose a normalized Gaussian with unit area underneath, amplitude determined fully by the FWHM
A = 1/np.sqrt((2*np.pi*var))

# Throw the points from the function using the utility function from the gauss class
points = gauss.draw_1dGauss(mean,var,N)

# Set num of bins

# Use Numpy to do the binning
#    hist is the freq hist in nbins, bin_edges are the left edge coords of these bins
hist, bin_edges = np.histogram(points,nbins) 


# Vars we'll need below
bin_centers = (bin_edges[:-1] + bin_edges[1:])/2 # Take avg to get x-coords of bin centers
paramvec = [A,mean,sigma]                        # Put the input params for the Gaussian into a vec

params_out, covar_matrix = curve_fit(gauss.gaussFun, bin_centers, hist, p0=paramvec)

# Make a vec of how many points we'll want to evaluate func at to make smooth curve, and the x domain
ncurvepoints = 1000
numsigmawidth = 5
x = np.linspace(mean-numsigmawidth*sigma,mean+numsigmawidth*sigma,ncurvepoints)

# Make the label for the fit Gaussian curve using output values
fit_gauss_label = plotutils.info_gaussian('Fit Gaussian',params_out[0],params_out[1],np.abs(params_out[2]))

# Plot the fit Gaussian curve (lw = linewidth)
fitcurve = plt.plot(x,gauss.gaussFun(x,*(params_out)),'r',lw=3,label=fit_gauss_label)

gaussLabel = " \n $\mu$ = " + str(mean) + ", A = " + str(np.round(A,3) ) + " , $\sigma$ = " + str(sigma)

# Make the title and x-y ax
es labels
plt.title('Gaussian with params:  '+gaussLabel ); plt.ylabel('Value'); plt.xlabel('x')

# Show the plot
plt.show()

sys.exit()
'''
