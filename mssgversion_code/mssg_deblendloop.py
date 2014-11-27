# 10-26-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import ipdb

### Declare what ellips to use for the objs
e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

e2ain=0.0
e2bin=0.0

### Whether to show the fits (setting this higher will show more things, current possibilities: 0, 1, 2)
pltflag = 0

############## Single shot
'''
outfile = "dumpfile"
e1ain = 0
e1bin = 0
print " ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin
# Make system call to fitting code
os.system('python mssg_deblendingTest.py  --e1a ' + str(e1ain) + ' --e2a ' + str(e2ain) + ' --e1b ' + str(e1bin) + ' --e2b ' + str(e2bin) + ' --outfile ' + outfile + ' --plotflag '+ str(pltflag) )
sys.exit()
'''


for numfiles in [0,1,2,3,4]:
### Run over ellips
    for e1bin in e1b_range:
        outfile = "deblendsOutput/deblResults_file"+str(numfiles)+"_e1b_in_" + str(e1bin)         
        with open(outfile+".txt", 'w') as f:
            f.write("{}\t".format('e1a_in'))
            f.write("{}\t".format('e1a_unbl'))
            f.write("{}\t".format('e1a_unblerr'))
            f.write( "{}\t".format('e1a_debl') )
            f.write( "{}\t".format('e1a_deblerr') )
            f.write( "{}\t".format('e1a_unbl - e1a_ein') )
            f.write( "{}\t".format('e1a_debl - e1a_ein') )
            
            f.write("{}\t\t".format('e1b_in'))
            f.write("{}\t".format('e1b_unbl'))
            f.write("{}\t".format('e1b_unblerr'))
            f.write( "{}\t".format('e1b_debl') )
            f.write( "{}\t".format('e1b_deblerr') )
            f.write( "{}\t".format('e1b_unbl - e1b_ein') )
            f.write( "{}\n".format('e1b_debl - e1b_ein') )

        for e1ain in e1a_range:

            print " ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin
        # Make system call to fitting code
            os.system('python mssg_deblendingTest.py  --e1a ' + str(e1ain) + ' --e2a ' + str(e2ain) + ' --e1b ' + str(e1bin) + ' --e2b ' + str(e2bin) + ' --outfile ' + outfile + ' --plotflag '+ str(pltflag))
