# 9-22-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import noiseLibrary
import  mssg_singleObjectNoiseStudy
import numpy as np

numtrials = 10

# SNR_range = [100,40,30,20,15,10,5]

SNR_range = [200,40,20]

e1in = 0.0
e2in = 0.0

for snr in SNR_range:
    print " *********** We're doing SNR = " , snr
    os.system('python mssg_singleObjectNoiseStudy.py  --niter ' +str(numtrials) + ' --snr ' + str(snr) + ' --e1 ' + str(e1in) + ' --e2 ' + str(e2in) )

######################################################################################################

#    Flux_range = [1e6,5e5,1e5,1e4,1e3,1e2]     # Full range
# Flux_range = [1e3,1e2]    

# texp = 6900 # seconds
# sbar = 26.8 # sky photons per second per pixel
# im = noiseLibrary.draw_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,

