# 10-24-2014
# MSSG


import matplotlib.pyplot as plt
import sys
import os
import noiseLibrary
import  mssg_singleObjectNoiseStudy
import numpy as np
import scipy.optimize as opt
 
e1bstr = str(0.5)


e1file = 'deblendsOutput/deblendingTests_e1b_in_' + e1bstr + '.txt'

f = open(e1file)
i = 0

vece1a_in = []
vece1a_exact = []
vece1a_debl = []
vece1a_resid = []

for line in f.readlines():
    splitline = line.strip().split("\t")
    if i > 0:
        vece1a_in.append(float(splitline[0]))
        vece1a_exact.append( float(splitline[1]))
        vece1a_debl.append( float(splitline[2]))
        vece1a_resid.append (float(splitline[3]))
    i += 1
    print i

print vece1a_in 
print vece1a_exact 
print vece1a_debl
print vece1a_resid 

f.close()

######## e1 
plt.figure()
xlimit = 0.4
#plt.scatter( vece1a_in, vece1a_exact   )
#plt.scatter( vece1a_in, vece1a_debl, color = 'g'   )
plt.scatter( vece1a_in, vece1a_resid, color = 'r'   )
#plt.errorbar( e1inRange, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
#plt.scatter( e1inRange, e1mean , color = 'r' )    
#xlimit += 0.1
#plt.xlim( -xlimit, xlimit)
#plt.plot([x1,x2],[y1,y2])    
plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
plt.savefig("resid_e1aBlended-e1aExact_vs_e1aIn_e1b_" + e1bstr + ".png")
plt.show()


#e1vec = np.genfromtxt(e1file)
#e2vec = np.genfromtxt(e2file)

#print e1vec 
