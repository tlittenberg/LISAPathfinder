import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from matplotlib import cm
import itertools

#3D
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import art3d

#2D Hist
import matplotlib.patches as patches
from matplotlib.path import Path

import matplotlib.colors
from matplotlib.colors import LogNorm
import copy
import colormap as colormap

#Confidence
import scipy as sp
import scipy.stats

import os, sys


H = 8.315000e-01            # Height of spacecraft [m]
xsc = np.zeros(8)           #initializing x array 
ysc = np.zeros(8)	    #initializing y array
xsc[0] =     -9.260000e-01 # x coordinate of spacecraft bottom deck corner 1 [m] 'SC_BOT_CORNER_1_X': 
ysc[0] =     -2.168000e-01 # y coordinate of spacecraft bottom deck corner 1 [m] 'SC_BOT_CORNER_1_Y': 
xsc[1] =     -9.260000e-01 # x coordinate of spacecraft bottom deck corner 2 [m] 'SC_BOT_CORNER_2_X': 
ysc[1] =      2.048000e-01 # y coordinate of spacecraft bottom deck corner 2 [m] 'SC_BOT_CORNER_2_Y': 
xsc[2] =     -5.263000e-01 # x coordinate of spacecraft bottom deck corner 3 [m] 'SC_BOT_CORNER_3_X': 
ysc[2] =     8.970000e-01  # y coordinate of spacecraft bottom deck corner 3 [m] 'SC_BOT_CORNER_3_Y': 
xsc[3] =     5.163000e-01  # x coordinate of spacecraft bottom deck corner 4 [m  'SC_BOT_CORNER_4_X': 
ysc[3] =     8.970000e-01  # y coordinate of spacecraft bottom deck corner 4 [m] 'SC_BOT_CORNER_4_Y': 
xsc[4] =     9.160000e-01  # x coordinate of spacecraft bottom deck corner 5 [m] 'SC_BOT_CORNER_5_X': 
ysc[4] =     2.048000e-01  # y coordinate of spacecraft bottom deck corner 5 [m] 'SC_BOT_CORNER_5_Y': 
xsc[5] =     9.160000e-01  # x coordinate of spacecraft bottom deck corner 6 [m] 'SC_BOT_CORNER_6_X': 
ysc[5] =     -2.168000e-01 # y coordinate of spacecraft bottom deck corner 6 [m] 'SC_BOT_CORNER_6_Y': 
xsc[6] =     5.163000e-01  # x coordinate of spacecraft bottom deck corner 7 [m] 'SC_BOT_CORNER_7_X': 
ysc[6] =     -9.090000e-01 # y coordinate of spacecraft bottom deck corner 7 [m] 'SC_BOT_CORNER_7_Y': 
xsc[7] =     -5.263000e-01 # x coordinate of spacecraft bottom deck corner 8 [m] 'SC_BOT_CORNER_8_X': 
ysc[7] =     -9.090000e-01 # y coordinate of spacecraft bottom deck corner 8 [m] 'SC_BOT_CORNER_8_Y': 


############# FUNCTIONS #############

### Funtion transforms around the 'origin', or the specific spacecraft coordinte, makes all y values equal
def transform(xs,ys,xorigin, yorigin, gotovector, index):
    
    xprime = []
    yprime = []
    
    for i in index: # Must index this way or the program will attempt to transform NaN
        r = np.sqrt((xs[i] - xorigin)**2 + (ys[i] - yorigin)**2) #finds distance away from sc (Space Craft) origin
        magvec = np.sqrt(gotovector[0]**2 + gotovector[1]**2)    #find the magnitude of the pointing vector
        #print(r)
        unitvec = [gotovector[0]/magvec,gotovector[1]/magvec]    #finds the unit vector 
        
        #normvec = np.linalg.norm(gotovector)
        xprime.append(r * unitvec[0]  + xorigin)                 #from the origin, adds x position
        yprime.append(r * unitvec[1]  + yorigin)		 #from the origin, adds y position
    return xprime, yprime
        

def findmb(x1,y1,x2,y2):
    m = (y2-y1)/(x2-x1)
    b = y1 - m*x1
    return m, b


def getfinite(df,i):
        facenumna = df.where(df['face'] == i)        #Makes places where face =! i NaN
        facenum = facenumna[np.isfinite(facenumna['impactnum'])] #Makes places where face =! dissapear 
        #print('len facenumber = %s'%(len(facenum['impactnum'])))
        index = np.asarray(list(facenum.index.values))                       #Lets you skip NaN stuff 
        if len(facenum['impactnum']) == 1:
            facenum = df[df.impactnum != 1]
            index = np.asarray(list(facenum.index.values))
            #facenum.drop(facenum.index[0,:])
            #print facenum
            #print('############ Dropped Len 1 !!!!!!!!!!!!!!')
            #facenum = facenum[np.isfinite(facenum['impactnum'])] #Makes places where face =! dissapear 
        return facenum, index
####### In the future, make N based on area #######

def patch3d_LPF(df,N, scale, mednmom,currun, dirpath, cmap):
    faces = np.arange(0,10) #Cycles through faces on the LPF
    alpha = 1 
    ec = 'white'#'auto' #edgecolor 
    lw = .02 
    #Parameterizing Visuals
    fig = plt.figure(figsize = (10,10))                #size of figure
    ax = fig.add_subplot(1,1,1, projection = '3d')     #add subplot with equal axes
    ax.set_aspect('equal')
   
    bot = False
    

    #Arbitrary Limit for finding Mistakes
    lim = 5 
    #ax.set_xlim(-lim,lim)
    #ax.set_ylim(-lim,lim)
    
    Ltotal = xsc[5] - xsc[0]                        #Total length of LPF
    Wtotal = ysc[2] - ysc[7]                        #Total Height of LPF
    
    ndensity = N/Ltotal                             #Bins / Unit Length
    binsheight = int(Wtotal * ndensity)             #Bins along the Y axis on top and bottom
    binsy = int(H*N/Ltotal)                         #Number of Bins in the y (z) direction    
    

    #X Y Z lines
    linex = liney = linez = np.arange(-1,2)
    zeros = np.zeros(len(linex))
    ax.plot(linex,zeros,zeros)
    ax.plot(zeros,liney,zeros)
    ax.plot(zeros,zeros,linez)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    

    indicies = []
    for i in faces:
        facenumna = df.where(df['face'] == i)                    #Makes places where face =! i NaN
        facenum = facenumna[np.isfinite(facenumna['impactnum'])] #Makes places where face =! dissapear
        indicies.append(len(list(facenum.index.values)))         #Lets you skip NaN stuff 
    #colormap stuff

    if 'log' in scale:   
        vmax = max(indicies)
        vmin = .5#1e-13
        norm = LogNorm(vmin = vmin, vmax = vmax)
        #filename = '%s/lpf_flat_log_%s.png'%(dirpath,mednmom)
        #filename = '%s/lpf_flat_log_%s.png'%(rundirpath,currun)
    else:
        vmin = 0
        vmax = max(indicies)#/N
        norm = matplotlib.colors.Normalize(vmin, vmax)
        #filename = '%s/lpf_flat_lin_%s.png'%(dirpath,mednmom)
        #filename = '%s/lpf_flat_lin_%s.png'%(rundirpath,currun)
    cmap = cmap
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap))
    my_cmap.set_bad(my_cmap(0))
   
    #Facecolors
    sidecolor = '#FF8C00'
    colortop = 'navy'
    
    for i in faces: #Loops through Faces to create flattened LPF
        #print('facenumber = %s, count = %s' %(i,count))
        #count += 1
        #print(i)
        facenum, index = getfinite(df,i)

        # This is the top, keep initial conditions, x = x, y = y ... 
        if i ==  9: 
            z = H   # Zposition            
	    #to Shift graphing Position, Must shift everything
            xpush = 0                      #Shift Parameter for x
            ypush = 0                      #Shift Parameter for y

            minfacex = xsc[0] + xpush
            maxfacex = xsc[5] + xpush
        
            minfacey = ysc[7] + ypush
            maxfacey = ysc[2] + ypush
         
	    #input data
            xs = facenum['xloc'] + xpush 
            ys = facenum['yloc'] + ypush  

            #Creates Histogram (NxN), Xedges (N), and Yedges (N)
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [N,binsheight],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            
	    #transforms the Histogram so it can be graphed
            Hist = Hist.T 
            
            mtright, btright = findmb(xsc[3],ysc[3],xsc[4], ysc[4])
            mbright, bbright = findmb(xsc[6],ysc[6],xsc[5], ysc[5])
            mbleft, bbleft = findmb(xsc[0],ysc[0],xsc[7], ysc[7])
            mtleft, btleft = findmb(xsc[1],ysc[1],xsc[2], ysc[2])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], yedges[t]  ,  z),
                        (xedges[i+1]   ,     yedges[t] , z),
                        (xedges[i+1],        yedges[t+1] , z),
                        (xedges[i],       yedges[t+1]   , z))]
                    if (xedges[i+1] > xsc[3]):
                        if not (mtright*xedges[i] + btright < yedges[t]) and not (mbright * xedges[i] + bbright > yedges[t+1]):
                            ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#

                    elif (xedges[i+1] < xsc[2]):
                        if not (mbleft*xedges[i+1] + bbleft > yedges[t+1]) and not (mtleft * xedges[i+1] + btleft < yedges[t]):
                            ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
                    else:            
                        ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
            
        #This is the bottom, flip initial conditions, x = x, y = -y, z = z 
        elif i ==  8: #Checked
            z = 0    #Z position
            
            facenum8 = facenum['impactnum']

            xpush = 0                      #Shift in the x direction
            ypush = 0 #H + (ysc[2] - ysc[7])  #Shift in the y direction

            minfacex = xsc[0] + xpush      
            maxfacex = xsc[5] + xpush
        
            maxfacey = -1 * ysc[7] + ypush
            minfacey = -1 * ysc[2] + ypush
        
            xbins = np.linspace(minfacex,maxfacex,N)
            ybins = np.linspace(minfacey,maxfacey,N)
            
            xs = facenum['xloc']
            ys = facenum['yloc']  	    
	    
            #Create Histogram
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [N,binsheight], 
                                         range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
        
            mtright, btright = findmb(xsc[3],ysc[3],xsc[4], ysc[4])
            mbright, bbright = findmb(xsc[6],ysc[6],xsc[5], ysc[5])
            mbleft, bbleft = findmb(xsc[0],ysc[0],xsc[7], ysc[7])
            mtleft, btleft = findmb(xsc[1],ysc[1],xsc[2], ysc[2])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], yedges[t]  ,  z),
                        (xedges[i+1]   ,     yedges[t] , z),
                        (xedges[i+1],        yedges[t+1] , z),
                        (xedges[i],       yedges[t+1]   , z))]
                    if (xedges[i+1] > xsc[3]):
                        if not (mtright*xedges[i] + btright < yedges[t]) and not (mbright * xedges[i] + bbright > yedges[t+1]):
                            ax.add_collection3d(Poly3DCollection(verts, alpha = alpha,edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#

                    elif (xedges[i+1] < xsc[2]):
                        if not (mbleft*xedges[i+1] + bbleft > yedges[t+1]) and not (mtleft * xedges[i+1] + btleft < yedges[t]):
                            ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
                    else:            
                        ax.add_collection3d(Poly3DCollection(verts, alpha = alpha,edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
        #This side is rotated 
        elif i == 3:
	    
	    #creates the vector pointing from vertex 4 to vertex 3, the base of side 3
            basevectorx = xsc[3] - xsc[4]
            basevectory = ysc[3] - ysc[4]
            basevector = [basevectorx,basevectory]
            
	    #Length of the Base
            width = np.sqrt(basevectorx**2+basevectory**2)
            binsx =  int(width * ndensity)                 # bins based on area
           
            #Bins are not exactly the same, but they are pretty close 
            lenbinsx = width/binsx
            lenbinsy = H/binsy
            #print('lenbins',lenbinsx,lenbinsy)

            #point that plot is turning around
            xorigin = xsc[4]          
            yorigin = ysc[4]    
             
            maxfacex = width #xsc[4]         
            minfacex = 0     #xsc[4] - width 
        
            minfacey = 0              
            maxfacey = H              
            
	    #vector points towards transformation
            gotovector = [1,0]
       	    
	    #Data to be Transformed 
            xin = facenum['xloc']
            yin = facenum['yloc'] 
           
            #Dummy Data
            #xin = (xsc[3])*np.ones(2000)
            #yin = ysc[3]*np.ones(2000)
            #index = np.arange(0,2000)

            #transforms data to y = yorigin 
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            #print(np.shape(xprime),np.shape(facenum['zloc'])) 
            
            #Plots the transformed Data
            #plt.scatter(xin,yin, color = 'red')
            #plt.scatter(xprime,yprime, color = 'green')
            
            xs = xprime - xorigin 
            ys = facenum['zloc']
            
            #Dummy Ys
            #ys = (H/3)*np.ones(2000) # facenum['zloc'] 
            
            #Creates Histogram in Easy (X,Z) reference frame
            Hist,xedges,yedges = np.histogram2d(xs, ys,bins = [binsx,binsy],
                    range = [[minfacex, maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
            #print(Hist) 

            #Flips Matrix Left to Right 
            Hist = np.fliplr(Hist)

            #print('Hist shape = %s, xedges = %s, yedges = %s'%(np.shape(Hist), np.shape(xedges), np.shape(yedges))) 
            #vector perpendicular to the base of the side 
            
            m, b = findmb(xsc[3], ysc[3], xsc[4],ysc[4])

            xedges = np.linspace(xsc[3], xsc[4],len(xedges), endpoint = True) 

            #print('Hist shape = %s, xedges = %s, yedges = %s'%(np.shape(Hist), np.shape(xedges), np.shape(yedges))) 
            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], m * xedges[i] + b, yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b, yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b, yedges[t+1]),
                        (xedges[i]      , m * xedges[i]     + b, yedges[t+1]))]
                    #print(Hist[t,i])
                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
             

        #this side is transformed, like 3
        elif i == 7: ## Checked ##
	    
            #creates the vector pointing from vertex 0 to vertex 7, the base of side 7
            basevectorx = xsc[7] - xsc[0]
            basevectory = ysc[7] - ysc[0]
            basevector = [basevectorx,basevectory]
            
	    #Length of the Base
            width = np.sqrt(basevectorx**2+basevectory**2)
            binsx =  int(width * ndensity)                 # bins based on area
           
            #Bins are not exactly the same, but they are pretty close 
            lenbinsx = width/binsx
            lenbinsy = H/binsy
            #print('lenbins',lenbinsx,lenbinsy)

            #point that plot is turning around
            xorigin = xsc[7]          
            yorigin = ysc[7]        
             
            maxfacex = width #xsc[4]         
            minfacex = 0     #xsc[4] - width 
        
            minfacey = 0              
            maxfacey = H              
            
	    #vector points towards transformation
            gotovector = [1,0]
       	    
	    #Data to be Transformed 
            xin = facenum['xloc']
            yin = facenum['yloc'] 
           
            #Dummy Data
            #xin = (xsc[7])*np.ones(20)
            #yin = ysc[7]*np.ones(20)
            #index = np.arange(0,20)

            #transforms data to y = yorigin 
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            
            xs = xprime - xorigin 
            ys = facenum['zloc']
            
            #Dummy Ys
            #ys = (0)*np.ones(20) # facenum['zloc'] 
            
            #Creates Histogram in Easy (X,Z) reference frame
            Hist,xedges,yedges = np.histogram2d(xs, ys,bins = [binsx,binsy],
                    range = [[minfacex, maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T

            #Flips Hist Left to right
            Hist = np.fliplr(Hist)
             
            xedges = np.linspace(xsc[0], xsc[7],len(xedges)) 
            m, b = findmb(xsc[0], ysc[0], xsc[7],ysc[7])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], m * xedges[i]     + b , yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b , yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b , yedges[t+1]),
                        (xedges[i]      , m * xedges[i]     + b , yedges[t+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
        
       #This side is transformed like 1
        elif i == 5: #Checked 
            
            binsx = int(width*ndensity) # bins based on area
            
            #base vector, pointing from 6 to 5
            basevectorx = xsc[5] - xsc[6]
            basevectory = ysc[5] - ysc[6]
            basevector = [basevectorx,basevectory]
            
            #width of base
            width = np.sqrt(basevectorx**2+basevectory**2)
       
            lenbinsx = width/binsx
            lenbinsy = H/binsy
            #print('lenbins',lenbinsx,lenbinsy)
            
            xpush = 0
            ypush = 0

            #i should probably change these to make them useful 
            minfacex = 0 #xsc[1] + xpush #xsc[1]  #+ width
            maxfacex = width  #xsc[1] + width + xpush        #xsc[1]  #width
        
            minfacey = 0      + ypush
            maxfacey = H      + ypush
             
            binsx = int(width * ndensity) #int((width*n / ltotal)) # bins based on area
            
            #print(binsx,binsy)
            
            #point that plot is turning around
            xorigin = xsc[6]
            yorigin = ysc[6]
            
            #direction transforming to, unit vector 
            gotovector = [1,0]

            #data to be transformed, currently dummy data
            xin = facenum['xloc']
            yin = facenum['yloc']

            #dummy data for tests
            #xin = (xsc[5])*np.ones(20) #+ xsc[2] + width    #facenum['xloc']
            #yin = [ysc[5]]*np.ones(20) #+ ysc[2]        #facenum['yloc']
            #index = np.arange(0,20)

            #transform data
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            
            #transformed data, figure out why xorigin must be added
            xs = xprime - xorigin 
            ys = facenum['zloc']

            #dummy data for y
            #ys = (H/3)*np.ones(20)#facenum['zloc']
            
            #create hist and edges from transformed data
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
            

            xedges = np.linspace(xsc[6], xsc[5],len(xedges)) 
            m, b = findmb(xsc[6], ysc[6], xsc[5],ysc[5])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], m * xedges[i] + b  , yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b   , yedges[t]),
                        (xedges[i+1]    , m * xedges[i+1]   + b   , yedges[t+1]),
                        (xedges[i]      , m * xedges[i]     + b   , yedges[t+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw , facecolor = my_cmap(norm(Hist[t,i]))))#
        
        
        ##this side is transformed like 5
        elif i == 1:

            #base vector, pointing from 1 to 2
            basevectorx = xsc[2] - xsc[1]
            basevectory = ysc[2] - ysc[1]
            basevector = [basevectorx,basevectory]
            
            #width of base
            width = np.sqrt(basevectorx**2+basevectory**2)
       
            lenbinsx = width/binsx
            lenbinsy = H/binsy
            #print('lenbins',lenbinsx,lenbinsy)
            xpush = 0
            ypush = 0

            #i should probably change these to make them useful 
            minfacex = 0 #xsc[1] + xpush #xsc[1]  #+ width
            maxfacex = width  #xsc[1] + width + xpush        #xsc[1]  #width
        
            minfacey = 0      + ypush
            maxfacey = H      + ypush
             
            binsx = int(width * ndensity) #int((width*n / ltotal)) # bins based on area
            
            #print(binsx,binsy)
            
            #point that plot is turning around
            xorigin = xsc[2]
            yorigin = ysc[2]
            
            #direction transforming to, unit vector 
            gotovector = [1,0]

            #data to be transformed, currently dummy data
            xin = facenum['xloc']
            yin = facenum['yloc']

            #dummy data for tests
            #xin = (xsc[1])*np.ones(20) #+ xsc[2] + width    #facenum['xloc']
            #yin = [ysc[1]]*np.ones(20) #+ ysc[2]        #facenum['yloc']
            #index = np.arange(0,20)

            #transform data
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            #plotting some stuff rn
            #plt.scatter(xin,yin,color = 'green')
            #plt.scatter(xprime,yprime, color = 'yellow')
            
            #transformed data, figure out why xorigin must be added
            xs = xprime - xorigin 
            ys = facenum['zloc']

            #dummy data for y
            #ys = (0)*np.ones(20)#facenum['zloc']
            
            #create hist and edges from transformed data
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T

            #Flips Histogram Left to Right
            Hist = np.fliplr(Hist) 
             
            xedges = np.linspace(xsc[1], xsc[2],len(xedges)) 
            m, b = findmb(xsc[1], ysc[1], xsc[2],ysc[2])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i] , m * xedges[i] + b  , yedges[t]),
                        (xedges[i+1], m * xedges[i+1]   + b   , yedges[t]),
                        (xedges[i+1], m * xedges[i+1]   + b   , yedges[t+1]),
                        (xedges[i]  , m * xedges[i]     + b   , yedges[t+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw , facecolor = my_cmap(norm(Hist[t,i]))))#
        
        ### Y is actually Z
        elif i == 2: # Checked 
            
            minfacex = xsc[2] 
            maxfacex = xsc[3]
        
            maxfacey = H 
            minfacey = 0 
            
            width = xsc[3] - xsc[2]

            binsx = int(ndensity * width) # bins based on area
           
 
            # Dummy Data 
            #xs = xsc[2]*np.ones(20) 
            #ys = H*np.ones(20)      
            
            xs = facenum['xloc'] 
            ys = facenum['zloc'] 

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            
            Hist = Hist.T 
            
            #print ('Hist = %s, xedges = %s, yedges = %s'%(np.shape(Hist), np.shape(xedges), np.shape(yedges)))
            
            m, b = findmb(xsc[2], ysc[2], xsc[3],ysc[3])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], m * xedges[i] + b, yedges[t]),
                        (xedges[i+1], m * xedges[i+1] + b, yedges[t]),
                        (xedges[i+1], m * xedges[i+1] + b, yedges[t+1]),
                        (xedges[i], m * xedges[i] +b , yedges[t+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw , facecolor = my_cmap(norm(Hist[t,i]))))#
            
        ### Y is actually Z
        elif i == 6: 
            
            minfacex = xsc[7] 
            maxfacex = xsc[6] 
        
            minfacey = 0 
            maxfacey = H 
            
            width = xsc[6] - xsc[7]

            binsx = int((width / Ltotal) * N) # bins based on area
           
            #Dummy
            #xs = xsc[6]*np.ones(20)
            #ys = H * np.ones(20)
            
            xs = facenum['xloc'] 
            ys = facenum['zloc'] 

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T 

            m, b = findmb(xsc[7], ysc[7], xsc[6],ysc[6])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((xedges[i], m * xedges[i] + b, yedges[t]),
                        (xedges[i+1], m * xedges[i+1] + b, yedges[t]),
                        (xedges[i+1], m * xedges[i+1] + b, yedges[t+1]),
                        (xedges[i], m * xedges[i] +b , yedges[t+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
       
       
        ### X is actually Z
        
        elif i == 4:                            # Checked, parameterized correctly
            z = xsc[5]
            xpush = 0
            ypush = 0 

            maxfacex = H 
            minfacex = 0 
           
            minfacey = ysc[0] 
            maxfacey = ysc[1] 
            
            width = maxfacey - minfacey

            binsx = int(ndensity * width) # bins based on area
            
            #xs = H*np.ones(20)
            #ys = ysc[4] * np.ones(20)

            xs = facenum['zloc'] 
            ys = facenum['yloc'] 

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsy,binsx], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
         
            m, b = findmb(ysc[5], xsc[5], ysc[4],xsc[4])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((m * yedges[t] + b, yedges[t]      , xedges[i]),
                        (m * yedges[t+1]     + b, yedges[t+1]    , xedges[i]),
                        (m * yedges[t+1]     + b, yedges[t+1]    , xedges[i+1]),
                        (m * yedges[t]       + b, yedges[t]      , xedges[i+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
        # X is actually Z 
        elif i == 0:               #Parameterized Correctly, Check done
            z = xsc[0]

            minfacex = 0 
            maxfacex = H 
           
            minfacey = ysc[0] 
            maxfacey = ysc[1] 
            
            width = maxfacey - minfacey

            binsx = int((width / Ltotal) * N) # bins based on area
            
            #xs = H * np.ones(20)
            #ys = ysc[0] * np.ones(20) 
            xs = facenum['zloc'] 
            ys = facenum['yloc'] 

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsy,binsx], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
        
            m, b = findmb(ysc[0], xsc[0], ysc[1],xsc[1])

            for t in range(len(yedges)-1):
                for i in range(len(xedges)-1):
                    verts = [((m * yedges[t] + b, yedges[t]      , xedges[i]),
                        (m * yedges[t+1]     + b, yedges[t+1]    , xedges[i]),
                        (m * yedges[t+1]     + b, yedges[t+1]    , xedges[i+1]),
                        (m * yedges[t]       + b, yedges[t]      , xedges[i+1]))]

                    ax.add_collection3d(Poly3DCollection(verts, alpha = alpha, edgecolor  = ec, linewidth = lw, facecolor = my_cmap(norm(Hist[t,i]))))#
        
    dirnamebot = '%s/%s_3dlpf_mom%s_run%s_bot'%(dirpath,scale,mednmom, currun)

    dirnametop = '%s/%s_3dlpf_mom%s_run%s_top'%(dirpath,scale,mednmom, currun)

    filestart = 'lpf_%s3d'%(scale)
        #filename = '%s/lpf_flat_log_%s.png'%(rundirpath,currun)
    #else:
        #filename = '%s/lpf_flat_lin_%s.png'%(dirpath,mednmom)
        #filename = '%s/lpf_flat_lin_%s.png'%(rundirpath,currun)
    #plt.tight_layout()
    #plt.savefig(filename)
    if not os.path.exists(dirnamebot):
        os.mkdir(dirnamebot)
    if not os.path.exists(dirnametop):
        os.mkdir(dirnametop)
    
    step = 15 
    #print(dirnametop)
    #print(dirnamebot)
    
    ax.set_title('3D Lisa Pathfinder Log Scale \n Time = %s s, Momentum = %s kg m/s'%(currun, mednmom), y=1.08)#, size = 'small')
    for ii in xrange(0,360,step):
        print('gif bottom angle = %s'%(ii))
        #elev += step 
        ax.view_init(elev = -15, azim = ii)
       
        plt.savefig("%s/%s_%s.png" %(dirnamebot,filestart,ii))
        #plt.clf()

    ax.set_title('3D Lisa Pathfinder Log Scale \n Time = %s s, Momentum = %s kg m/s'%(currun, mednmom), y = 1.08)#, size = 'small')
    for ii in xrange(0,360+step,step):
        print('gif top  angle = %s'%(ii))
        #elev += step 
        ax.view_init(elev = 15, azim = ii)
        plt.savefig("%s/%s_%s.png" %(dirnametop,filestart,ii))
        #plt.clf()
    plt.close('all')
    #plt.show()
            


