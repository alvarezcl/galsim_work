# Start: 10-24-2014
# MSSG

# Takes the files on disk and plots them up 

import matplotlib.pyplot as plt
import sys
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--subdir", default="deblendsOutput/", help="output text filename")
args = parser.parse_args()

############ Read in file
fname = args.subdir + 'deblendingTests_5runs.txt'
fitdat = np.loadtxt(fname)

fnumvec = fitdat[:,0]

# Initze all a vecs
e1a_in = fitdat[:,1]
e1a_unbl = fitdat[:,3] #;      e1a_unblerr = []
e1a_debl = fitdat[:,4] # ;       # e1a_deblerr = []
e1a_unblresid =  e1a_unbl - e1a_in 
e1a_deblresid =  e1a_debl - e1a_in 

e1b_in = fitdat[:,5]
e1b_unbl = fitdat[:,7] #;      e1b_unblerr = []
e1b_debl = fitdat[:,8] # ;       # e1b_deblerr = []
e1b_unblresid =  e1b_unbl - e1b_in 
e1b_deblresid =  e1b_debl - e1b_in 

fnum =  fitdat[:,0]

print 'fnum = ', fnum 
print 'e1a_in = ' ,e1a_in 
print 'e1a_unbl  = ' ,e1a_unbl 
print 'e1a_debl = ', e1a_debl
print 'e1a_unblresid = ',e1a_unblresid

print 'e1b_in = ' ,e1b_in 
print 'e1b_unbl  = ' ,e1b_unbl 
print 'e1b_debl = ', e1b_debl
print 'e1b_unblresid = ',e1b_unblresid

numfiles = 5

#    for i in xrange(0,numfiles):

e1a_range = [0.5, 0.25, 0, -0.25, -0.5]
e1b_range = [0.5, 0.25, 0, -0.25, -0.5]

for e1bin in e1b_range:
    e1bstr = str(e1bin)

    vece1a_in = []  
    vece1a_unbl = []
    vece1a_debl = []
    vece1a_unblresid = []  
    vece1a_deblresid = []  
    vece1a_unblerr = []  
    vece1a_deblerr = []  

    for e1ain in e1a_range:
        i = np.where(np.logical_and(e1a_in == e1ain, e1b_in == e1bin))
        print 'e1bin, e1ain,  i = ', e1bin, e1ain, i 
#        print 'fitdat[i, 3] = ', fitdat[i, 3] 

        vece1a_in.append( e1a_in[i].mean() )
        vece1a_unbl.append( e1a_unbl[i].mean() )
        vece1a_debl.append( e1a_debl[i].mean() )
        vece1a_unblresid.append (e1a_unblresid[i].mean()  )
        vece1a_deblresid.append (e1a_deblresid[i].mean() )
        vece1a_unblerr.append( e1a_unbl[i].std()  )
        vece1a_deblerr.append( e1a_debl[i].std() ) 

######## e1 a
    print 'vece1a_unbl = ', vece1a_unbl
    plt.figure()
    xlimit = 0.4
    plt.scatter( e1a_range, vece1a_unbl, color = 'g'   )
#    plt.scatter( e1a_in[i], vece1a_debl, color = 'g'   )
#plt.scatter( vece1a_in, vece1a_resid, color = 'r'   )

   
#    gline = plt.errorbar( vece1a_in, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1a" )
#    bline = plt.errorbar( vece1a_in, vece1a_deblresid, vece1a_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1a" ) 
   
    plt.title("Resids for fits with e1bin = " +e1bstr )
    plt.legend() # (handles=[gline,bline])
    plt.xlabel('$e_{1a in}$',fontsize=18)
    plt.ylabel('$e_{1a fit}-e_{1a in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
    plt.show()

#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")

sys.exit()

e1a_in = 0
e1a_unbl = 0 ;      e1a_unblerr = 0
e1a_debl = 0 ;       e1a_deblerr = 0
e1a_unblresid = 0 ; e1a_deblresid = 0

#  e1a_unbl += vece1a[where vece1b_unbl== e1bin]





numf = 5
for filenum in xrange(0,5):
    prependname = args.subdir + 'deblResults_file_'+str(filenum)+'_e1b_in_'

    for e1bin in e1b_range:
        e1bstr = str(e1bin)

        e1file = prependname + e1bstr + '.txt'
        print e1file

        f = open(e1file)
        print "Plotting up ",e1file

# Read in the whole file, and put each column entry into a vec elt
        i = 0
        for line in f.readlines():
            print "\n\n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> i = ", i
            
            splitline = line.strip().split("\t")

  #      print splitline
#        print 'splitline[7] = ', splitline[7]
 #       print 'splitline[8] =',splitline[8]
#        if i>1:
 #           sys.exit()

            if i > 0:  # First line is column headers
                e1a_in +=float(splitline[0])
                e1a_unbl +=float(splitline[1])
                e1a_debl +=float(splitline[3])
                e1a_unblresid +=float(splitline[5])
                e1a_deblresid +=float(splitline[6])
            
                vece1b_in.append(float(splitline[7]))
            # Skip 8 because there were 2 tabs in this file
                vece1b_unbl.append( float(splitline[9]))
                vece1b_unblerr.append( float(splitline[10]))
                vece1b_debl.append( float(splitline[11]))
                vece1b_deblerr.append( float(splitline[12]))
                vece1b_unblresid.append (float(splitline[13]))
                vece1b_deblresid.append (float(splitline[14]))
                i += 1

    vece1a_in.append(e1a_in / numf )
    vece1a_unbl.append(e1a_unbl / numf )
    vece1a_debl.append( e1a_debl / numf )
    vece1a_unblresid.append (e1a_unblresid / numf  )
    vece1a_deblresid.append (e1a_deblresid / numf )
    vece1a_unblerr.append( float(splitline[2]))
    vece1a_deblerr.append( float(splitline[4]))
                
    vece1b_in.append(float(splitline[7]))
            # Skip 8 because there were 2 tabs in this file
    vece1b_unbl.append( float(splitline[9]))
    vece1b_unblerr.append( float(splitline[10]))
    vece1b_debl.append( float(splitline[11]))
    vece1b_deblerr.append( float(splitline[12]))
    vece1b_unblresid.append (float(splitline[13]))
    vece1b_deblresid.append (float(splitline[14]))
    i += 1
                
#print vece1a_in 
#print vece1a_unbl 
#print vece1a_debl

    f.close()

######## e1 a
    plt.figure()
    xlimit = 0.4
#plt.scatter( vece1a_in, vece1a_unbl   )
#plt.scatter( vece1a_in, vece1a_debl, color = 'g'   )
#plt.scatter( vece1a_in, vece1a_resid, color = 'r'   )
    gline = plt.errorbar( vece1a_in, vece1a_unblresid, vece1a_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1a" )
    bline = plt.errorbar( vece1a_in, vece1a_deblresid, vece1a_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1a" ) 
#plt.errorbar( e1inRange, e1mean , yerr= [e/sqrtnumt for e in e1sigma] ,ecolor='g',linestyle=' ' )
#plt.scatter( e1inRange, e1mean , color = 'r' )    
#xlimit += 0.1
#plt.xlim( -xlimit, xlimit)
#plt.plot([x1,x2],[y1,y2])    
    plt.title("Resids for fits on: " + e1file)
    plt.legend() # (handles=[gline,bline])
    plt.xlabel('$e_{1a in}$',fontsize=18)
    plt.ylabel('$e_{1a fit}-e_{1a in}$',fontsize=18)
    plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")
    plt.show()

######## e1b 
print vece1b_in 
print vece1b_unbl 
print vece1b_debl


plt.figure()
xlimit = 0.4
gline = plt.errorbar( vece1b_in, vece1b_unblresid, vece1b_unblerr,  ecolor='g',linestyle=' ', label = "Resid for unblended fit for e1b" )
bline = plt.errorbar( vece1b_in, vece1b_deblresid, vece1b_deblerr,  ecolor='b',linestyle=' ', label = "Resid for deblended fit for e1b" ) 
plt.title("Resids for fits on e1b " )
plt.legend() # (handles=[gline,bline])
plt.xlabel('$e_{1b in}$',fontsize=18)
plt.ylabel('$e_{1b fit}-e_{1b in}$',fontsize=18)
plt.axhline( 0,color='k',linestyle='-', linewidth= 2)     
#    plt.savefig("resid_e1aBlended-e1aUnbl_vs_e1aIn_e1b_" + e1bstr + ".png")
plt.show()
