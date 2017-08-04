import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from matplotlib import cm
import itertools
import matplotlib.lines as mlines

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

from decimal import Decimal
import os, sys

##Parametrizing Spacecraft Geometry, changed to array below because this was annoying to type
Geo = {'SC_H' : 8.315000e-01, # Height of spacecraft [m]
        'SC_BOT_CORNER_1_X': -9.260000e-01, # x coordinate of spacecraft bottom deck corner 1 [m]
        'SC_BOT_CORNER_1_Y': -2.168000e-01, # y coordinate of spacecraft bottom deck corner 1 [m]
        'SC_BOT_CORNER_2_X': -9.260000e-01, # x coordinate of spacecraft bottom deck corner 2 [m]
        'SC_BOT_CORNER_2_Y':  2.048000e-01, # y coordinate of spacecraft bottom deck corner 2 [m]
        'SC_BOT_CORNER_3_X': -5.263000e-01, # x coordinate of spacecraft bottom deck corner 3 [m]
        'SC_BOT_CORNER_3_Y': 8.970000e-01,  # y coordinate of spacecraft bottom deck corner 3 [m]
        'SC_BOT_CORNER_4_X': 5.163000e-01,  # x coordinate of spacecraft bottom deck corner 4 [m
        'SC_BOT_CORNER_4_Y': 8.970000e-01,  # y coordinate of spacecraft bottom deck corner 4 [m]
        'SC_BOT_CORNER_5_X': 9.160000e-01,  # x coordinate of spacecraft bottom deck corner 5 [m]
        'SC_BOT_CORNER_5_Y': 2.048000e-01,  # y coordinate of spacecraft bottom deck corner 5 [m]
        'SC_BOT_CORNER_6_X': 9.160000e-01,  # x coordinate of spacecraft bottom deck corner 6 [m]
        'SC_BOT_CORNER_6_Y': -2.168000e-01, # y coordinate of spacecraft bottom deck corner 6 [m]
        'SC_BOT_CORNER_7_X': 5.163000e-01,  # x coordinate of spacecraft bottom deck corner 7 [m]
        'SC_BOT_CORNER_7_Y': -9.090000e-01, # y coordinate of spacecraft bottom deck corner 7 [m]
        'SC_BOT_CORNER_8_X': -5.263000e-01, # x coordinate of spacecraft bottom deck corner 8 [m]
        'SC_BOT_CORNER_8_Y': -9.090000e-01, # y coordinate of spacecraft bottom deck corner 8 [m]
        'EOM_RB_X' : 1.303672e-03,  # X of S/C CoM in M Frame => this defines the B frame [m]
        'EOM_RB_Y' : 2.370536e-03,  # Y of S/C CoM in M Frame => this defines the B frame [m]
        'EOM_RB_Z' : 4.913037e-01   # Z of S/C CoM in M Frame => this defines the B frame [m]
        }

H = 8.315000e-01            #Height of spacecraft [m]
xsc = np.zeros(8)           #initializing x array 
ysc = np.zeros(8)           #initializing y array
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

## Import Data ##
home = False

############# FUNCTIONS #############

### Funtion transforms around the 'origin', or the specific spacecraft coordinte, makes all y values equal, used for Flat LPF
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
        yprime.append(r * unitvec[1]  + yorigin)                 #from the origin, adds y position
    return xprime, yprime
        

# Rotates Matrix, Used for Flat LPF 
def DoRotation(xspan, yspan, RotRad): #(xspan -> flattened x coords, yspan -> flattend y coords)
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],  #Definition of rotation matrix
                          [-np.sin(RotRad), np.cos(RotRad)]])

    x, y = np.meshgrid(xspan, yspan)
    #plt.scatter(x,y,color = 'orange', alpha = .2)  
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

#Returns only finite values of the dataframe, gets rid of NaN from masking, Only works with newer numpy!
def getfinite(df,i):
        facenumna = df.where(df['face'] == i)        #Makes places where face =! i NaN
        facenum = facenumna.dropna(axis=0, how='any')#Makes places where NaN dissapear 
        
	#print('len facenumber = %s'%(len(facenum['impactnum']))) 
	if len(facenum['impactnum']) == 1: 		#if there is only one value in a face, just drop it. need [] or len > 1 for hist
            facenum = df[df.impactnum != 1] 		#drops all value where impact num =1 , impact num always = 1 ;) 
            index = np.asarray(list(facenum.index.values))
            
	    #print(facenum) #Should be empty. Just making sure
            #print('############ Dropped Len 1 !!!!!!!!!!!!!!')

        else:
	    index = np.asarray(list(facenum.index.values)) #get index 
	return facenum, index


#### Makes the Flattened LPF, the thing that you fold up ####
# N for each side based on area given an initial N #

def flatten_LPF(df, N, scale, mednmom,currun, cmap, dirpath, rundirpath):
    faces = np.arange(0,10) 			     #Cycles through faces on the LPF
    alpha = 0 					     #Facecolor of patches, currently transparent

    #Parameterizing Visuals
    fig = plt.figure(figsize = (10,10))              #size of figure
    ax = fig.add_subplot(1,1,1, aspect = 'equal')    #add subplot with equal axes
    ax.set_xlim(-2,2)                                #xlim
    ax.set_ylim(-2,4)                                #ylim
    

    #Arbitrary Limit for finding Mistakes
    #lim = 5 
    #ax.set_xlim(-lim,lim)
    #ax.set_ylim(-lim,lim)
    
    Ltotal = xsc[5] - xsc[0]                        #Total length of LPF
    Wtotal = ysc[2] - ysc[7]                        #Total Width of LPF (across top)
    
    #These names are confusing I'm sorry :// 
    ndensity = N/Ltotal                             #Bins per unit length, Bins / Unit Length
    binsheight = int(Wtotal * ndensity)             #Bins along the Y axis on top and bottom
    binsy = int(H*N/Ltotal)                         #Number of Bins in the y (z) direction    
    
    #Used to find mistakes 
    #Scattering an origin/lines
    #xline = np.arange(-1,1,.1)
    #yline = np.arange(-1,1,.1)
    #plt.plot(xline,np.zeros(len(xline)),color = 'black')
    #plt.plot(np.zeros(len(yline)),yline, color = 'black')
    
    #Old way of Parameterizing Spacecraft, I am too lazy to change
    points = [1,2,3,4,5,6,7,8,1]                    
    bottom = []
    top = []
    xy = []
    for p in points:   
        Xs = Geo['SC_BOT_CORNER_%s_X'%(p)]
        Ys = Geo['SC_BOT_CORNER_%s_Y'%(p)]
        xy.append([Xs,Ys])
        top.append((Xs,Ys,H))
        bottom.append((Xs,Ys,0))
     
    #Find Length of sides for normalizing 
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
        filename = '%s/lpf_flat_log_%s_%s.png'%(dirpath,currun,mednmom)
    #    filenamerun = '%s/lpf_flat_log_%s.png'%(rundirpath,currun)
    else:
        vmin = 0
        vmax = max(indicies)/N
        norm = matplotlib.colors.Normalize(vmin, vmax)
        filename = '%s/lpf_flat_lin_%s_%s.png'%(dirpath,currun,mednmom)
    #    filenamerun = '%s/lpf_flat_lin_%s.png'%(rundirpath,currun)
    cmap = cmap
    #Sets bad values (ie: log(0) to the lowest value on the map)
    my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap))
    my_cmap.set_bad(my_cmap(0))
   
    #Facecolors
    sidecolor = '#FF8C00' #orange
    colortop = 'navy'     #navy

    count = 0
    lennums = []

    #Honestly this function is just so ugly that I don't want to comment it .... 
    for i in faces: #Loops through Faces to create flattened LPF
        #print('facenumber = %s, count = %s' %(i,count))
        count += 1
            

        facenum, index = getfinite(df,i)
	lennums.append((len(facenum['impactnum']))*100/len(df['impactnum']))
	
	#X and Z switched
        if i == 0:               #Parameterized Correctly, Check done
	    #The left most face
            z = xsc[0]
            xpush = xsc[0]-H #Places the X correctly on the flattened LPF, does not change data
            ypush = 0

            minfacex = 0 + xpush
            maxfacex = H + xpush
           
            minfacey = ysc[0] + ypush
            maxfacey = ysc[1] + ypush
            
            width = maxfacey - minfacey

            binsx = int((width / Ltotal) * N) # bins based on area
            binsy = int(H*ndensity)
	    xs = facenum['zloc'] + xpush
            ys = facenum['yloc'] + ypush

	    xs = np.asarray(xs)
	    ys = np.asarray(ys) 
	    #print(np.shape(xs), np.shape(ys), np.shape([binsx,binsy]))
	 
	    Hist,xedges,yedges = np.histogram2d(xs,ys, bins = [binsy,binsx], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T 
            
	    #Makes Patch
	    xyside0 = [[minfacex,minfacey], [maxfacex, minfacey],[maxfacex, maxfacey], [minfacex, maxfacey]]
            path0 = Path(xyside0)
            patch0 = patches.PathPatch(path0, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch0)
	    #Plots Color
            ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patch0, clip_on=True)
        

        ##this side is transformed like 5
        elif i == 1:

            #base vector, pointing from 1 to 2
            basevectorx = xsc[2] - xsc[1]
            basevectory = ysc[2] - ysc[1]
            basevector = [basevectorx,basevectory]
            
            #width of base
            width = np.sqrt(basevectorx**2+basevectory**2)
       
            xpush = 0
            ypush = 0

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

            #transform data, flattens the side so that there are no angles
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            #plotting some stuff rn
            #plt.scatter(xin,yin,color = 'green')
            #plt.scatter(xprime,yprime, color = 'yellow')
            
            #transformed data, figure out why xorigin must be added
            xs = xprime - xorigin 
            ys = facenum['zloc']

            #dummy data for y
            #ys = (0)*np.ones(20)#facenum['zloc']
            np.asarray(xs)
	    np.asarray(ys) 
            
	    #create hist and edges from transformed data
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
            
            #plot unrotated hist, where is it?? 
            #ux,uy = np.meshgrid(xedges,yedges) 
            #ax.pcolormesh(ux,uy,hist, cmap = 'summer')
             
            #find angles between sides
            perpbase = [-basevector[1],basevector[0]]  #vector perpendicular to the base
            
            vec1 = basevector
            vec2 = [-1,0]
           
            theta = 1*np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))) #angle between sides
            
            #print(np.degrees(theta)) 

            farvec_constant =  (H/(np.linalg.norm(basevector)))
            farvec = np.multiply(farvec_constant, perpbase)
            
            #creating vertecies for patch
            xmax = xsc[1]  +farvec[0]
            ymax = ysc[1]  +farvec[1]
            cornx = xsc[2] +farvec[0]
            corny = ysc[2] +farvec[1]
            xyside1 = [[xsc[1],ysc[1]], [xsc[2], ysc[2]],[cornx,corny],[xmax,ymax]]
             
            #places patch in the right spot
            xplace = xorigin - H * np.sin(theta) 
            yplace = yorigin - H * np.cos(theta) 
            
            #patch stuff
            path1 = Path(xyside1)
            patch1 = patches.PathPatch(path1, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch1)
            
            #rotate hist, rotate the histogram so that the sides are angled
            x,y = DoRotation(xedges,yedges,(theta))
            ax.pcolormesh(x+xplace,y+yplace,Hist, norm = norm, cmap = my_cmap, alpha = 1, clip_path = patch1, clip_on = True)
            
            #Plot some other stuff 
            #plt.scatter(x[0],y[0], color = 'red', alpha = .2)
            #plt.scatter(x[1],y[1],color = 'green', alpha = .2)
            #print(x[0][0],y[0][0])
	
        ### Y is actually Z
        elif i == 2:  #Topmost side that is not the bottom (8)
            xpush = 0
            ypush = 0 #ysc[2] + H
            
            minfacex = xsc[7] + xpush
            maxfacex = xsc[6] + xpush
        
            maxfacey = H  + ypush
            minfacey = 0  + ypush 
            
            width = xsc[6] - xsc[7]

            binsx = int((width / Ltotal) * N) # bins based on area
           

            xs = facenum['xloc'] + xpush
            ys = facenum['zloc'] + ypush

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T 
            
	    #Flips histogram up and down (visual purposes)
	    Hist = np.flipud(Hist)
            xyside2 = [[minfacex,ysc[2]], [maxfacex, ysc[2]],[maxfacex, ysc[2]+H], [minfacex, ysc[2]+H]]

            path2 = Path(xyside2)
            patch2 = patches.PathPatch(path2, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch2)
            
            xedges  = np.linspace(xsc[2], xsc[3], len(xedges))
            yedges = np.linspace(ysc[2], ysc[2] + H, len(yedges))
            
            #Plots the hist
            ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patch2, clip_on=True)
	
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
            #xin = (xsc[3])*np.ones(20)
            #yin = ysc[3]*np.ones(20)
            #index = np.arange(0,20)

            #transforms data to y = yorigin 
            xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
            
            #print(np.shape(xprime),np.shape(facenum['zloc'])) 
            
            #Plots the transformed Data
            #plt.scatter(xin,yin, color = 'red')
            #plt.scatter(xprime,yprime, color = 'green')
            
            xs = xprime - xorigin 
            ys = facenum['zloc']
            
            #Dummy Ys
            #ys = (H/2)*np.ones(20) # facenum['zloc'] 
            np.asarray(xs)
	    np.asarray(ys) 
            #Creates Histogram in Easy (X,Z) reference frame
            Hist,xedges,yedges = np.histogram2d(xs, ys,bins = [binsx,binsy],
                    range = [[minfacex, maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
         
	    #Plotting stuff to find errors 
            #fig3 = plt.figure(figsize = (10,10))                #size of figure
            #ax3 = fig3.add_subplot(1,1,1, aspect = 'equal')#, projection = '3d')    #add subplot with equal axes
            #Plots untransformed Hist
            #ux,uy = np.meshgrid(xedges,yedges)
            #ax3.pcolormesh(ux,uy,Hist,norm = norm, cmap = cmap)
            #ax3.scatter(xin,yin,facenum['zloc'])
            #ax3.scatter(xs,ys)


            #vector perpendicular to the base of the side 
            perpbase = [basevector[1],-basevector[0]]
           
            #Find angle between vectors 
            vec1 = basevector
            vec2 = [1,0]
            
            #Angle between vectors, radians
            theta = 1*np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))) 
            #print(np.degrees(theta)) 

            farvec_constant = (H/(np.linalg.norm(basevector))) 
            farvec = np.multiply(farvec_constant, perpbase)     #Unit vector point towards top corner
            
            xmax = xsc[3]+farvec[0]                             #X position of top right
            ymax = ysc[3]+farvec[1]                             #Y position of top right
            cornx = xsc[4] +farvec[0]                           #X position of bot right
            corny = ysc[4] +farvec[1]                           #Y position of bot right
            xyside3 = [[cornx,corny],[xsc[4],ysc[4]], [xsc[3], ysc[3]],[xmax,ymax],[cornx,corny]] #corners for patch
            
            offsetx = .062009  #these were constants that kept coming out, I don't know why
            offsety = -.0873899

	    #Trig to figure out placement on flattened LPF
            xplace = xsc[4] + H * np.sin(theta)  #ysc[4]+ysc[3]  
            yplace = ysc[4] - H * np.cos(theta) #ysc[4] + ysc[4]#- H * np.sin(theta)  #xsc[4]+xsc[3] 

            path3 = Path(xyside3)
            patch3 = patches.PathPatch(path3, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch3)
            
            #Rotates Matrix by theta radians
            x,y = DoRotation(xedges,yedges,(-theta))
            ax.pcolormesh(x+xplace,y+yplace,Hist,norm = norm, cmap = my_cmap, clip_path = patch3, clip_on = True)
            
        ### X is actually Z
        
        elif i == 4:                            # Checked, parameterized correctly
            z = xsc[5]
            xpush = xsc[5] + H
            ypush = 0 #ysc[4]

            maxfacex = 0 + xpush
            minfacex = -1*H + xpush
           
            minfacey = ysc[0] + ypush
            maxfacey = ysc[1] + ypush
            
            width = maxfacey - minfacey

            binsx = int((width / Ltotal) * N) # bins based on area

            xs = -1 * facenum['zloc'] + xpush
            ys = facenum['yloc'] + ypush

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsy,binsx], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
        
            xyside4 = [[minfacex,minfacey], [maxfacex, minfacey],[maxfacex, maxfacey], [minfacex, maxfacey]]

            path4 = Path(xyside4)
            patch4 = patches.PathPatch(path4, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch4)

            ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patch4, clip_on=True)

        #This side is transformed like 1
        elif i == 5:
 
            binsx = int(width*ndensity) # bins based on area
            
            #base vector, pointing from 6 to 5
            basevectorx = xsc[5] - xsc[6]
            basevectory = ysc[5] - ysc[6]
            basevector = [basevectorx,basevectory]
            
            #width of base
            width = np.sqrt(basevectorx**2+basevectory**2)
       
            xpush = 0
            ypush = 0
	    
	    #Pretend that this side is not rotated
            minfacex = 0 
            maxfacex = width  
        
            minfacey = 0      + ypush
            maxfacey = H      + ypush
             
            binsx = int(width * ndensity) # bins based on area
            
            #point that plot is turning around
            xorigin = xsc[6]
            yorigin = ysc[6]
            
            #direction transforming to, unit vector 
            gotovector = [1,0]

            #data to be transformed, currently dummy data
            xin = facenum['xloc']
            yin = facenum['yloc']

            #dummy data for tests
            #xin = (xsc[6])*np.ones(20) #+ xsc[2] + width    #facenum['xloc']
            #yin = [ysc[6]]*np.ones(20) #+ ysc[2]        #facenum['yloc']
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
            np.asarray(xs)
	    np.asarray(ys) 
            #create hist and edges from transformed data
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
            
            #plot unrotated hist, where is it?? 
            #ux,uy = np.meshgrid(xedges,yedges) 
            #ax.pcolormesh(ux,uy,Hist, cmap = 'summer')
 
            #find angles between sides
            perpbase = [-basevector[1],basevector[0]]  #vector perpendicular to the base
            
            vec1 = basevector
            vec2 = [-1,0]
           
            theta = 1*np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))) #angle between sides
            
            #print(np.degrees(theta)) 

            farvec_constant =  (H/(np.linalg.norm(basevector)))
            farvec = np.multiply(farvec_constant, perpbase)
            
            #creating vertecies for patch
            xmax  = xsc[6]  - farvec[0]
            ymax  = ysc[6]  - farvec[1]
            cornx = xsc[5]  - farvec[0]
            corny = ysc[5]  - farvec[1]
            xyside5 = [[xsc[6],ysc[6]], [xsc[5], ysc[5]],[cornx,corny],[xmax,ymax], [xsc[6],ysc[6]]]
             
            #places patch in the right spot
            xplace = xorigin + H * np.sin(theta) 
            yplace = yorigin + H * np.cos(theta) 
            
            #patch stuff
            path5  = Path(xyside5)
            patch5 = patches.PathPatch(path5, facecolor=sidecolor, lw=2, alpha = alpha)
            ax.add_patch(patch5)
            
            #rotate hist
            x,y = DoRotation(-xedges,-yedges,(theta))
            ax.pcolormesh(x+xplace,y+yplace,Hist, norm = norm,cmap = my_cmap, alpha = 1, clip_path = patch5, clip_on = True)
            
            #Plot some other stuff 
            #plt.scatter(x[0],y[0], color = 'red', alpha = .2)
            #plt.scatter(x[1],y[1],color = 'green', alpha = .2)
            #print(x[0][0],y[0][0])
        
            #print(x[0][0],y[0][0])
        ### Y is actually Z
        elif i == 6: 
            xpush = 0
            ypush = 0 #ysc[6] - H  #ysc[6] + H
            
            minfacex = xsc[7] + xpush
            maxfacex = xsc[6] + xpush
        
            minfacey = 0 + ypush
            maxfacey = H + ypush 
            
            width = xsc[6] - xsc[7]

            binsx = int((width / Ltotal) * N) # bins based on area
           

            xs = facenum['xloc'] + xpush
            ys = facenum['zloc'] + ypush

            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [binsx,binsy], range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T 
             
            xyside6 = [[minfacex,ysc[7]-H], [maxfacex, ysc[7]-H],[maxfacex, ysc[7]], [minfacex, ysc[7]]]

            path6 = Path(xyside6)
            patch6 = patches.PathPatch(path6, facecolor=sidecolor, lw=2, alpha = alpha)
            
            ax.add_patch(patch6)
            
            xedges  = np.linspace(xsc[7], xsc[6], len(xedges))
            yedges = np.linspace(ysc[7]-H, ysc[7], len(yedges))
            
            #Plots the hist
            ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patch6, clip_on=True)

        #this side is transformed, like 3
        elif i == 7: 
            
            #creates the vector pointing from vertex 0 to vertex 7, the base of side 7
            basevectorx = xsc[0] - xsc[7]
            basevectory = ysc[0] - ysc[7]
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
            
            #print(np.shape(xprime),np.shape(facenum['zloc'])) 
            
            #Plots the transformed Data
            #plt.scatter(xin,yin, color = 'red')
            #plt.scatter(xprime,yprime, color = 'green')
            
            xs = xprime - xorigin 
            ys = facenum['zloc']
            
            #Dummy Ys
            #ys = (0)*np.ones(20) # facenum['zloc'] 
            np.asarray(xs)
	    np.asarray(ys) 
            
	    #Creates Histogram in Easy (X,Z) reference frame
            Hist,xedges,yedges = np.histogram2d(xs, ys,bins = [binsx,binsy],
                    range = [[minfacex, maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
           
            #Plots untransformed Hist
            #ux,uy = np.meshgrid(xedges,yedges)
            #ax.pcolormesh(ux,uy,Hist,cmap = 'summer')

            #vector perpendicular to the base of the side 
            perpbase = [basevector[1],-basevector[0]]
           
            #Find angle between vectors 
            vec1 = basevector
            vec2 = [1,0]
            
            #Angle between vectors, radians
            theta = 1*np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))) 
            #print(np.degrees(theta)) 

            farvec_constant = (H/(np.linalg.norm(basevector))) 
            farvec = np.multiply(farvec_constant, perpbase)     #Unit vector point towards top corner
            
            xmax =  xsc[0] -farvec[0]                             #X position of top right
            ymax =  ysc[0] -farvec[1]                             #Y position of top right
            cornx = xsc[7] -farvec[0]                           #X position of bot right
            corny = ysc[7] -farvec[1]                           #Y position of bot right
            xyside7 = [[cornx,corny],[xsc[7],ysc[7]], [xsc[0], ysc[0]],[xmax,ymax],[cornx,corny]] #corners for patch

            xplace = xsc[7] - H * np.sin(theta)  #ysc[4]+ysc[3]  
            yplace = ysc[7] + H * np.cos(theta) #ysc[4] + ysc[4]#- H * np.sin(theta)  #xsc[4]+xsc[3] 

            path7 = Path(xyside7)
            patch7 = patches.PathPatch(path7, facecolor = sidecolor , lw=2, alpha = alpha)
            ax.add_patch(patch7)
            
            #Rotates Matrix by theta radians
            x,y = DoRotation(xedges,-yedges,(-theta))
            ax.pcolormesh(x+xplace,y+yplace,Hist,norm = norm, cmap = my_cmap, clip_path = patch7, clip_on = True)
            
        #This is the bottom, flip initial conditions, x = x, y = -y, z = z 
        elif i ==  8: 
            z = 0    #Z position

            xpush = 0                      #Shift in the x direction
            ypush = ysc[2] + H - ysc[7]  #Shift in the y direction

            minfacex = xsc[0] + xpush      
            maxfacex = xsc[5] + xpush
        
            maxfacey = -1 * ysc[7] + ypush
            minfacey = -1 * ysc[2] + ypush
        
            xbins = np.linspace(minfacex,maxfacex,N)
            ybins = np.linspace(minfacey,maxfacey,N)
            
            xs = facenum['xloc']    + xpush
            ys = -1*facenum['yloc'] + ypush  #Flipped y because the bottom is viewed upside down
            xs = np.asarray(xs)
	    ys = np.asarray(ys) 
            #Create Histogram
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [N,binsheight], 
                                         range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            Hist = Hist.T
            
            #Creates Patch for bottom 
            xybot = []
            for i in range(len(xsc)):
                xybot.append([xsc[i]+xpush,ysc[i]+ypush])
            pathbot = Path(xybot)
            patchbot = patches.PathPatch(pathbot, facecolor=sidecolor, lw=2, alpha = alpha)
            patchbot1 = patchbot
            ax.add_patch(patchbot) 

            #Plots the hist, and gets cropped by the octogon
            ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patchbot, clip_on=True)
	
	
	# This is the top, keep initial conditions, x = x, y = y ... 
        elif i ==  9: 
            
	    z = H   # Zposition           
            
            #To Shift graphing Position, Must shift everything
            xpush = 0                      #Shift Parameter for x
            ypush = 0                      #Shift Parameter for y

            minfacex = xsc[0] + xpush
            maxfacex = xsc[5] + xpush
        
            minfacey = ysc[7] + ypush
            maxfacey = ysc[2] + ypush
         
            #Input data
            xs = facenum['xloc'] + xpush 
            ys = facenum['yloc'] + ypush  
	    
	    xs = np.asarray(xs)
	    ys = np.asarray(ys)
            #Creates Histogram (NxN), Xedges (N), and Yedges (N)
            Hist,xedges,yedges = np.histogram2d(xs,ys,bins = [N,binsheight],range = [[minfacex,maxfacex],[minfacey,maxfacey]])
            
            #Transforms the Histogram so it can be graphed
            Hist = Hist.T 
            
            #Creates the Octogon Patch 
            xytop = []
            for i in range(len(xsc)):
                xytop.append([xsc[i]+xpush,ysc[i]+ypush])
            pathtop = Path(xytop)
            patchtop = patches.PathPatch(pathtop, facecolor=colortop, lw=2, alpha = alpha)
            patchtop1 = patchtop
            ax.add_patch(patchtop) 

            #Plots the hist, and gets cropped by the octogon
            plottop =ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
                    clip_path = patchtop, clip_on=True)
        
            plt.colorbar(plottop) #Makes the colorbar for all graphs, normalization is the same 
    
    #Labels the facenumbers, 2 and 9 are missing because they are ugly when put on 
    marking = my_cmap(0)
    ax.annotate('0', xy=(-1.85, 0), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('1', xy=(-1.5, 1), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('3', xy=(1.5, 1), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('4', xy=(1.85, 0), color = marking)#, xytext=(xsc[0] - H -.5,0)) 
    ax.annotate('5', xy=(1.5, -1), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('6', xy=(0, -1.9), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('7', xy=(-1.6,-1), color = marking)#, xytext=(xsc[0] - H -.5,0))
    ax.annotate('8', xy=(0, 3.6), color = marking)#, xytext=(xsc[0] - H -.5,0))
    #ax.annotate('9', xy=(1.5, 1), color = marking)#, xytext=(xsc[0] - H -.5,0))
    
    #Makes Legend for percent hit of side numbers
    textstr = 'Percent Hit \n'
    for i in range(10):
	textstr += 'Face %i = %i'%(i, lennums[i])
	if i != 9:
		textstr +='\n'
    props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

    ax.text(0.05, 0.97, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props) 
    ax.set_title('LPF Impact Location %s\n Median Momentum = %s kg m/s, Time = %s s'%(scale, mednmom, currun), size = 'small')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    plt.close()

    #print('closed fig %s'%(filename))
            

            
### Just makes a 3D model of LPF, nothing scientific ###
createlpf = False
if createlpf:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    alpha = .75
    xypoints = []
    for i in range(7):
        xypoints.append((xsc[i],ysc[i]))

    points = [1,2,3,4,5,6,7,8,1]
    bottom = []
    top = []
    xy = []
    for p in points:   
        Xs = Geo['SC_BOT_CORNER_%s_X'%(p)]
        Ys = Geo['SC_BOT_CORNER_%s_Y'%(p)]
        xy.append([Xs,Ys])
        top.append((Xs,Ys,H))
        bottom.append((Xs,Ys,0))
    
    sideslist = []
    for p in range(8):
        sides = (top[p],top[p+1],bottom[p+1],bottom[p])
        #print(sides)
        vertsside = [sides]
        
        sideslist.append([top[p],top[p+1],bottom[p+1],bottom[p]])
        #facecolors = [cm.jet(Histtop)]
        sideart = ax.add_collection3d(Poly3DCollection(vertsside, alpha = .75))#, facecolor = colormap))
    #print(ax.get_children()) 
    vertstop = [top]
    vertsbot = [bottom]
    topart = ax.add_collection3d(Poly3DCollection(vertstop, alpha = .75))
    botart = ax.add_collection3d(Poly3DCollection(vertsbot, alpha = .75))


    linex = liney = linez = np.arange(-1,2)
    zeros = np.zeros(len(linex))
    ax.plot(linex,zeros,zeros)
    ax.plot(zeros,liney,zeros)
    ax.plot(zeros,zeros,linez)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


### Finds confidence interval (ignore the mean)
def mean_confidence_interval(data, confidence):
 
    #First sorts data in order
    data = np.sort(data, axis = None)
    N = len(data)

    condown =((1-confidence)/2)                 #percentage down
    conup   = (confidence + ((1-confidence)/2)) #percentage up
    
    nup = int(N * conup)    			#index value of confidence up
    ndown = int(N * condown) 			#Index value of confidence down
    
    confidence_up = data[nup]
    confidence_down = data[ndown]
    
    return confidence_down, confidence_up

#Plots NOISE not error, 90% confidence value of NOISE
def ploterror(mednmom,currun, dirpath, rundirpath, homedir):    
    freqs   = []
    
    Snx     = []
    Sny     = []
    Snz     = []
    Sntheta = []
    Sneta   = []
    Snphi   = []


    filenums = np.arange(50,100) #Skips first 50 files, Run in Reasons
    
    for i in filenums:
            impactname = homedir + '/runs/run_e_%s/psd.dat.%s'%(currun,i)#'~/website/runs/run_e_%s/psd.dat.%s'%(currun,i)
            dfnoise = pd.read_csv(impactname, header = None,delimiter = '\s+',
                names = ['freq','Snx','Sny','Snz','Sntheta','Sneta','Snphi'],skiprows = 1, skipfooter = 1, engine = 'python')
            freqs.append(dfnoise['freq'].values)
            Snx.append(dfnoise['Snx'].values)
            Sny.append(dfnoise['Sny'].values)
            Snz.append(dfnoise['Snz'].values)
            Sntheta.append(dfnoise['Sntheta'].values)
            Sneta.append(dfnoise['Sneta'].values)
            Snphi.append(dfnoise['Snphi'].values)

    freq = dfnoise['freq']
    freqsnow = freqs
    
    params = [Snx,Sny,Snz,Sntheta,Sneta,Snphi]
    names =  ['Snx', 'Sny', 'Snz', 'Sntheta', 'Sneta', 'Snphi']
    
    fig,axs = plt.subplots(figsize = (6,12),nrows = len(params), ncols = 2)#, sharex = True)
    t = 0
    errorcolor = '#F5BE46'
    linecolor  = '#E06959'
    fig.suptitle('PSD Noise \n Median Momentum = %s, Time = %s'%(mednmom,currun)) 
    for p in params:
            p = np.asarray(p)
            #p_mean = []
            p_errorup = []
            p_errordown = []
            for i in range(len(freq)):
                    #Grabs the first row in all files to get the first freq, should be 50 values
                    perrordown, perrorup= mean_confidence_interval(p[:,i], confidence=0.90)
                    p_errordown.append(perrordown), p_errorup.append(perrorup)
            
	    #t tells you what param, snx,sny ...
	    modelname = homedir + '/runs/run_e_%s/FD_model_%s.dat'%(currun,t) 
	    dfmodel = pd.read_csv(modelname, header = None,delimiter = '\s+',
                names = ['freq','noise','dontcare','dontknow'],skiprows = 0, skipfooter = 0, engine = 'python')
            
	    p_errordown = np.asarray(p_errordown)
            p_errorup = np.asarray(p_errorup)
            #print(p_errordown)
            ax0 = axs[t,0]
            ax0.set_ylabel(names[t])
	    #Scaling Axes 
            ax0.set_yscale('log')
            ax0.set_xscale('linear')
    	    ax0.set_xlim(min(freq), max(freq))
	    ax0.plot(dfmodel['freq'], dfmodel['noise'],'-', alpha = .75, color = 'teal', label = 'Model Noise')
            ax0.plot(freq, p_errorup/2,'-', color = linecolor, label = '90% Confidence Noise')
            ax0.plot(freq,p_errordown/2,color = linecolor)
            ax0.fill_between(freq, p_errordown/2, p_errorup/2, color = errorcolor)  
	    if t == 0: 
                ax0.set_title('Semilog Y', size = 'small')

            ax1 = axs[t,1]
            ax1.set_yscale('log')
	    ax1.set_xscale('log')
    	    ax1.set_xlim(min(freq), max(freq))
	    ax1.plot(dfmodel['freq'], dfmodel['noise'],'-', alpha = .75, color = 'teal', label = 'Model Noise')
            ax1.plot(freq, p_errorup/2,'-', color = linecolor, label = '90% Confidence Noise')
            ax1.plot(freq,p_errordown/2,color = linecolor)
            ax1.fill_between(freq, p_errordown/2, p_errorup/2, color = errorcolor) 
	    if t == 0: 
                ax1.set_title('Log Log', size = 'small')
	    
            t += 1 
    #Makes Legend
    modelline = mlines.Line2D([],[],color = 'teal', label = 'Model Noise') #model Error
    bounds = mlines.Line2D([],[],color=linecolor, label = '90% Confidence Noise') #90% confidence
    
    ax = axs[0,0]
    ax.legend(loc = 'lower left', fontsize = 'small')#handles = [modelline,bounds], labels = ['Model Error', '90% Confidence'],loc = 'lower left')#handles = [modelline,bounds],labels = ['Model Error', '90% Confidence'], loc='lower left')#(boundline,errorline),('90% Confidence', 'Model Noise'), 'upper right')
    ax = axs[5,0]
    ax.set_xlabel('Freq Hz')
    ax = axs[5,1]
    ax.set_xlabel('Freq Hz', size = 'small')
    plt.tight_layout(rect=[0,0.03,1,0.95])
    plt.savefig('%s/dofsnoise_currun%s_mednmom%s.png'%(dirpath,currun,mednmom)) 
    plt.close()

