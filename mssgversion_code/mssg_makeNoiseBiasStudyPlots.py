# 9-22-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import noiseLibrary
import  mssg_singleObjectNoiseStudy
import numpy as np
import scipy.optimize as opt

def linefunc(x, b, m):
    return m*x + b

SNR_range = [200,40,20]

########################## Making the plots

path = 'noisebiasoutput/'

e1val = 0.0
e2val = 0.0

e1inRange = [-0.1, -0.05, 0.0, 0.07, 0.2]
e2inRange = [-0.1, -0.05, 0.0, 0.07, 0.2]

# e1inRange = [-0.1]


os.system('rm '+path+'m1VsSNR.txt')
os.system('rm '+path+'m2VsSNR.txt')

# Loop over SNR values
for snr in SNR_range:
    print " snr  = ", snr
    e1mean = [] ; e2mean = []
    e1sigma = [] ; e2sigma = []


    # Loop over input values of e1
    for e1val in e1inRange:

        print " *********** We're doing SNR = " , snr

        # Get the filenames
        e1file = path + "noise_bias_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e1out.txt" 
        e2file = path + "noise_bias_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e2out.txt" 

        # Make data vector from the file
        e1vec = np.genfromtxt(e1file)
        e2vec = np.genfromtxt(e2file)
        
        e1mean.append(np.mean(e1vec) - e1val)
        e2mean.append(np.mean(e2vec) - e2val)
        
        e1sigma.append(np.std(e1vec) )
        e2sigma.append(np.std(e2vec) )
        
        print " e1mean = ", e1mean 
        print " e2mean = ", e2mean 
        
        print " e1sigma = ", e1sigma 
        print " e2sigma = ", e2sigma 

        
############# e1

    x0    = np.array([0.0, 0.0])
    results1 =  opt.curve_fit(linefunc, e1inRange, e1mean, x0)
    
    print results1
    
    b1 = results1[0][0]
    m1 = results1[0][1]
    m1err = np.sqrt(results1[1][1][1])

    print " slope1 = ", m1
    print " y-intercept1 = ", b1
    print " m1err^2  = ", m1err*m1err

############# e2

    x0    = np.array([0.0, 0.0])
    results2 =  opt.curve_fit(linefunc, e2inRange, e2mean, x0)
    
    print results2
    
    b2 = results2[0][0]
    m2 = results2[0][1]
    m2err = np.sqrt(results2[1][1][1])
    print " slope2 = ", m2
    print " y-intercept2 = ", b2
    print " m2err^2  = ", m2err*m2err

# Write to files
    with open(path+'m1VsSNR.txt', 'a') as f:
        f.write("{}\n".format(snr))
        f.write("{}\n".format(m1))
        f.write("{}\n".format(m1err))

    with open(path+'m2VsSNR.txt', 'a') as f:
        f.write("{}\n".format(snr))
        f.write("{}\n".format(m2))
        f.write("{}\n".format(m2err))


# Now plot these    
#    plt.figure()
    numtrials = 100
    sqrtnumt = np.sqrt(numtrials)
   
######## e1 
    plt.figure()
    xlimit = 0.4
    x1 = -xlimit ; y1 = m1*x1 + b1;
    x2 =  xlimit ; y2 = m1*x2 + b1;       
    plt.errorbar( e1inRange, e1mean , yerr= e1sigma ,ecolor='b',linestyle=' ' )
    plt.errorbar( e1inRange, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
    plt.scatter( e1inRange, e1mean , color = 'r' )    
    xlimit += 0.1
    plt.xlim( -xlimit, xlimit)
    plt.plot([x1,x2],[y1,y2])    
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
#    plt.show()
    plt.savefig("e1out_vs_SNR.png")

########### e2
    plt.figure()
    xlimit = 0.4
    x1 = -xlimit ; y1 = m2*x1 + b2;
    x2 =  xlimit ; y2 = m2*x2 + b2;
    plt.errorbar( e2inRange, e2mean , yerr= e2sigma ,ecolor='b',linestyle=' ' )
    plt.errorbar( e2inRange, e2mean , yerr= [e/sqrtnumt for e in e2sigma] ,ecolor='g',linestyle=' ' )
    plt.scatter( e2inRange, e2mean , color = 'r' )    
    xlimit += 0.1
    plt.xlim( -xlimit, xlimit)
    plt.plot([x1,x2],[y1,y2])    
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)    
#    plt.show()    
    plt.savefig("e2out_vs_SNR.png")
    
sys.exit()

################################################################3





# leftover: true_e1 = 0.0;    plt.axhline( true_e1,color='k',linestyle='-')


#plt.errorbar( SNR_range, e1mean , yerr= e1sigma ,ecolor='b',linestyle=' ' )
#plt.errorbar( SNR_range, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
#plt.scatter( SNR_range, e1mean , color = 'r' )
    
