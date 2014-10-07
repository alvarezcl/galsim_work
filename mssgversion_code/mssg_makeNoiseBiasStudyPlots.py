# 9-22-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import noiseLibrary
import  mssg_singleObjectNoiseStudy
import numpy as np
import scipy.optimize as opt

def func(x, b, m):
    return b + m*x 

SNR_range = [100,40,30,20,15,10,5]

SNR_range = [200]

e1mean = []
e2mean = []

e1sigma = []
e2sigma = []

########################## Making the plots

path = 'noisebiasoutput/'

e1val = 0.0
e2val = 0.0

e1inRange = [-0.1, -0.05, 0.0, 0.07, 0.2]

# e1inRange = [-0.1]

# Loop over input values of e1
for e1val in e1inRange:

# Loop over SNR values
    for snr in SNR_range:
        print " *********** We're doing SNR = " , snr

        # Get the filenames
        e1file = path + "noise_bias_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e1out.txt" 
        e2file = path + "noise_bias_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e2out.txt" 

        # Make data vector from the file
        e1vec = np.genfromtxt(e1file)
        e2vec = np.genfromtxt(e2file)
        
        e1mean.append(np.mean(e1vec))
        e2mean.append(np.mean(e2vec))
        
        e1sigma.append(np.std(e1vec))
        e2sigma.append(np.std(e2vec))
        
        print " e1mean = ", e1mean 
        print " e2mean = ", e2mean 
        
        print " e1sigma = ", e1sigma 
        print " e2sigma = ", e2sigma 

        
############# e1

x0    = np.array([0.0, 0.0])
print opt.curve_fit(func, e1inRange, e1mean, x0)


plt.figure()
numtrials = 100
sqrtnumt = np.sqrt(numtrials)

plt.errorbar( e1inRange, e1mean , yerr= e1sigma ,ecolor='b',linestyle=' ' )
plt.errorbar( e1inRange, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
plt.scatter( e1inRange, e1mean , color = 'r' )


#plt.errorbar( SNR_range, e1mean , yerr= e1sigma ,ecolor='b',linestyle=' ' )
#plt.errorbar( SNR_range, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
#plt.scatter( SNR_range, e1mean , color = 'r' )

true_e1 = 0.0

plt.axhline( true_e1,color='k',linestyle='-')
plt.axhline( 0,color='k',linestyle='-', linewidth= 2)



plt.show()

plt.savefig("e1out_vs_SNR.png")

sys.exit()

########### e2
plt.figure()

plt.errorbar( SNR_range, e2mean , yerr= e2sigma ,ecolor='b',linestyle=' ' )
plt.errorbar( SNR_range, e2mean , yerr= [e/sqrtnumt for e in e2sigma] ,ecolor='g',linestyle=' ' )
plt.scatter( SNR_range, e2mean , color = 'r' )

true_e2 = 0.0

plt.axhline( true_e2,color='k',linestyle='-')
plt.axhline( 0,color='k',linestyle='-', linewidth= 2)

plt.show()

plt.savefig("e2out_vs_SNR.png")
