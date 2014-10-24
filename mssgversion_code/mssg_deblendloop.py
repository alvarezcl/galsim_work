# 10-26-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import noiseLibrary
import numpy as np
import ipdb

e1a_range = [0.5, 0.25, 0.00, -0.25, -0.5]
#e2a_range = [0.5, 0.01, -0.5]

e1b_range = [0.5, 0, -0.5]

# e1bin = e1b_range[0]


e2ain=0
e2bin=0

### Run over ellips
for e1bin in e1b_range:
    outfile = "deblendsOutput/deblendingTests_e1b_in_" + str(e1bin)        
    with open(outfile+".txt", 'w') as f:
        f.write("{}\t".format('e1a_in'))
        f.write("{}\t".format('e1a_exact'))
        f.write( "{}\t".format('e1a_debl') )
        f.write( "{}\n".format('e1a_debl - e1a_exact') )

    for e1ain in e1a_range:

        print " ************************************************ We're doing e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin
        # Call fitting code
        os.system('python mssg_deblendingTest.py  --e1a ' + str(e1ain) + ' --e2a ' + str(e2ain) + ' --e1b ' + str(e1bin) + ' --e2b ' + str(e2bin) + ' --outfile ' + outfile)
