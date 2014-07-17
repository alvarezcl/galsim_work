# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 16:19:06 2014

@author: luis
"""

## Library of functions with useful formatting and plotting
## functions. 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from matplotlib import cm

# Return a string (with Latex formatting) of particular information
# in the following format: Type,Amplitude,mean,sigma
def info_gaussian(type,A,mu,sigma):    
    return type + ': $A=%.2f$, $\mu=%.2f$, $\sigma=%.2f$' % (A,mu,sigma)

# Return similar string with primary gaussian
def info_gaussian_Acp(type,Ar,mu,sigma):
    return type + ': $A_{cp}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$' % (Ar,mu,sigma)

# Return similar string with secondary gaussian
def info_gaussian_Acs(type,Ar,mu,sigma):
    return type + ': $A_{cs}=%.2f$, $\mu=%.2f$, $\sigma=%.2f$' % (Ar,mu,sigma)


# Plot 3d surface given any set of values, X, Y, Z    
def plot_3d(X,Y,Z):
    fig1 = plt.figure()
    ax1 = Axes3D(fig1)
    surf = ax1.plot_surface(X,Y,Z,cmap=cm.coolwarm)
    fig1.colorbar(surf,shrink=0.5,aspect=5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('f(x,y)')
    return fig1
    
# Plot contours of the functional values in Z    
def plot_contour(Z):    
    fig2 = plt.figure(2)
    im = pl.imshow(Z)
    cset = pl.contour(Z)
    pl.clabel(cset,inline=True)
    pl.colorbar(im)
    plt.title('Contour')
    return fig2

# Given a data array of x and y such that they correspond
# to a number of points, this draws a 3d Bar Graph corresponding
# to density of points.
def hist_3dBar_points(x,y,bins):
    fig = plt.figure()
    ax = Axes3D(fig)
    hist, xedges, yedges = np.histogram2d(x,y,bins)
    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.3*np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
    return fig
    
def hist_3dBar_binned(X,Y,H):
    fig = plt.figure()
    ax = Axes3D(fig)
    xpos = X.flatten()
    ypos = Y.flatten()
    zpos = np.zeros(X.shape[0]*Y.shape[1])
    dx = 0.3*np.ones_like(zpos)
    dy = dx.copy()
    dz = H.flatten()
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
    return fig
    

# Produce 2D histogram projection of the density of points corresponding
# to x-y.   
def hist_2dPlane(x,y,bins):
    H,xedges,yedges = np.histogram2d(x,y,bins,normed=False)
    X,Y = np.meshgrid(xedges,yedges)
    img = plt.imshow(H,interpolation ='none')
    #plt.grid(True)
    pl.colorbar(img)
    plt.show()
    return X,Y,H,img