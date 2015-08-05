# MSSG
# 8/4/2015

# Simple code to plot up outputs from deblendloopOverSep.py

from pylab import *

datawi = genfromtxt('deblendsOutput/sepstep_10StepsAcrossPixel.withinterp.txt')
datafi1 = genfromtxt('deblendsOutput/sepstep_10StepsAcrossPixel.withforceinterp.txt')
datafi2 = genfromtxt('deblendsOutput/sepstep_10StepsAcrossPixel.withforceinterp2.txt')
datafi3 = genfromtxt('deblendsOutput/sepstep_10StepsAcrossPixel.withforceinterp3.txt')
datan = genfromtxt('deblendsOutput/sepstep_10StepsAcrossPixel.nointerp.txt')

plot(datan[:,17], datan[:,10], "o")
plot(datawi[:,17], datawi[:,10], "o")
plot(datafi1[:,17], datafi1[:,10], "o")
plot(datafi2[:,17], datafi2[:,10], "o")
plot(datafi3[:,17], datafi3[:,10], "o")

show()
