import matplotlib.pyplot as plt
import numpy as np

#################################################### Function to draw gal figs
def makeplot(pltname, pltcontent, centre=(0,0) , colorbarlimit = 0, location = 'lower', mrkrsize = 10 ):


          plt.title(pltname)
          print " >>>>>> Plotting ", pltname          
          dummyarr = np.array([1, 2, 3])
          if type(pltcontent) != type(dummyarr):
              print pltcontent, " is not a numpy array, will convert"
              pltcontent = pltcontent.array 
          pltimg = plt.imshow( pltcontent , origin=location, interpolation='none' );
          if colorbarlimit==1:
                    pltimg.set_clim(-1,1)
          plt.colorbar()
          if (centre[0] !=0 and centre[1] !=0):
                    plt.plot(centre[0][0], centre[0][1],  color='g', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
                    plt.plot(centre[1][0], centre[1][1],  color='y', marker = 'o', markersize=mrkrsize, markeredgewidth=1)
                
          plt.show()
    
