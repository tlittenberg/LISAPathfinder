# microTools.py - a python module for doing analysis on LPF micrometeoroid data
"""
microTools is a python module for performing analysis of LPF micrometeoroite impact events
that have been identified and characterized using an MCMC tool. It contains functions to
read the chain files, produce plots, conduct population analysis, etc.

Ira Thorpe
2018-05-12
"""

# function to read chain files as output by MCMC tool
def readRawChain(chainDir,grs=1, burnIn=0.5, outDir='data'):
    """
    Function to read micrometeoroite MCMC chain files. Assumes directory structure 
    as in the Aug2017 runs directory on tsankawi. Output is a python pickle file 
    containing a dictionary with the relevant chain information. 
    Arguments
       chainDir = directory corresponding to the MCMC output is found
       grs = index of the grs for this chain (1 or 2)
       burnIn = fraction of chain to throw out for burn in
       outDir = output directory for pickle file (name will be generated automatically)
       
    Ira Thorpe
    2018-05-12 
    """
    
    # modules
    import numpy as np
    import os
    import pickle
    import string
    
    # find directory and get gps time
    base = os.path.basename(chainDir)
    gpsTime = float(base[len(base)-10:len(base)])
    
    # load impactChain
    impFile = chainDir +'/impactchain.dat'
    dat = np.loadtxt(impFile)
    N = np.shape(dat)[0]
    trim = int(float(N)*burnIn);
    dat = np.delete(dat, slice(0, trim), axis=0)
    
    # build into a dictionary
    data = {
    	'gps' : gpsTime,
    	'N' : np.shape(dat)[0],
    	't0' : dat[:,2], 
    	'Ptot' : dat[:,3], 
    	'lat' : 90-(np.arccos(dat[:,7])*180/np.pi), 
    	'lon' : np.mod(dat[:,8]*180/np.pi+180,360)-180,
    	'rx' : dat[:,10],
    	'ry' : dat[:,11],
    	'rz' : dat[:,12],
    	'face' : dat[:,9]}
    	
    # load log likelihood chain
    logLfile = chainDir +'/logLchain.dat'
    dat = np.loadtxt(logLfile)
    N = np.shape(dat)[0]
    trim = int(float(N)*burnIn);
    dat = np.delete(dat, slice(0, trim), axis=0)
    
    # compute detection fraction
    dfrac = np.sum(dat[:,0])/(np.shape(dat)[0])
    data['dfrac'] = dfrac
    
    # save data in processed directory
    pickle.dump(data,open(str(os.cwd)+'/' + outDir + '/' + str(int(gpsTime))+'_grs' + str(int(grs)) + '.pickle','wb'))

    # return data
    return data
    
# function to get spacecraft quaternions    
def getSCquats(gps,doText=False):
    """
    function to read SC quaternion file. Can either read a python binary file (faster, 
    default) or an ASCII text file (slower)
    
    Ira Thorpe
    2018-05-12
    """
    #import libraries
    import numpy as np, quaternion
    import os
    import pathlib
    
    # get current working directory
    p = pathlib.PurePath(os.getcwd())
    baseDir = str(p.parent)
    
    # load Quaternion data
    if doText:
        quatFile = baseDir +'/rawData/allQuats.txt'
        dat = np.loadtxt(quatFile)
    else :
        quatFile = baseDir + '/data/quats.npy'
        dat = np.load(quatFile)
        
    # separate out gps time (1st column) from quaternions (columns 2-5)
    allGPS = np.array(dat[...,0])
    allQuats = quaternion.as_quat_array(np.array(dat[...,[4,1,2,3]]))
    
    # find nearest gps time
    idxmin = (np.abs(allGPS-gps)).argmin()
    impQuat = allQuats[idxmin]

    # return the quaternion
    return impQuat
    

# function to locate impact and estimate area using healpix binning.
def findSkyAngles(data, CI=0.68, nside=32):
    """
    function to determine impact sky area using HEALPIX binning. Returns 1 sigma sky area in 
    square degrees and central point latitude and longitude. If dictionary passed to the 
    function has sun-frame angles in addition to SC-frame angles, it will operate on both.
    Arguments
        data = dictionary containing chain data
        CI = confidence interval for sky area
        nside = HEALPIX number of sides
    
    Ira Thorpe
    2018-05-24
    """
    #import libraries
    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt

    # Build the HEALPIX map
    npix = hp.nside2npix(nside)
    mp = np.arange(npix)

    # Convert data to HEALPIX
    dat_hp = hp.pixelfunc.ang2pix(nside, data['lon'], data['lat'], nest=False, lonlat=True)

    # Make the histogram
    bin_edges = np.arange(-0.5,npix+0.5,1.0)
    bin_centers = np.arange(0,npix,1.0)
    cnt_hp, bins = np.histogram(dat_hp,bin_edges)
    
    # Measure centroid and sky area
    cdf = np.cumsum(cnt_hp.astype('float'))/float(data['N'])        
    ilb = (np.abs(cdf-((1.0-CI)/2.0))).argmin()
    iub = (np.abs(cdf-(1.0-((1.0-CI)/2.0)))).argmin()
    imed = (np.abs(cdf-0.5)).argmin()
    area = 41253.0*float(iub-ilb)/float(npix)
    lon_c, lat_c = hp.pixelfunc.pix2ang(nside, imed, nest=False, lonlat=True)
    lon_c = np.mod(180+lon_c,360)-180
    
    # put back into data dictionary
    data['lat_c'] = lat_c
    data['lon_c'] = lon_c
    data['skyArea'] = area
    data['healPix'] = cnt_hp/float(data['N'])
    
    # if Sun angles are present, repeat for them
    if 'lon_sun' in data :
        # Convert data to HEALPIX
        dat_hp = hp.pixelfunc.ang2pix(nside, data['lon_sun'], data['lat_sun'], nest=False, lonlat=True)

        # Make the histogram
        cnt_hp, bins = np.histogram(dat_hp,bin_edges)
    
        # Measure sky area
        cdf = np.cumsum(cnt_hp.astype('float'))/float(data['N'])        
        ilb = (np.abs(cdf-((1.0-CI)/2.0))).argmin()
        iub = (np.abs(cdf-(1.0-((1.0-CI)/2.0)))).argmin()
        imed = (np.abs(cdf-0.5)).argmin()
        area = 41253.0*float(iub-ilb)/float(npix)
        lon_c, lat_c = hp.pixelfunc.pix2ang(nside, imed, nest=False, lonlat=True)
        lon_c = np.mod(180+lon_c,360)-180
        
        # put into dictionary
        data['lat_c_sun'] = lat_c
        data['lon_c_sun'] = lon_c
        data['skyArea_sun'] = area
        data['healPix_sun'] = cnt_hp/float(data['N'])

    
    # return dictionary
    return data

# function to convert angles from SC frame to Sun-center frame (in degrees)
def SCtoSun(data):
	"""
	funciton to convert angles from SC frame to sun-centered frame used by micrometeoroid 
	population models. 
	
	Ira Thorpe
	2018-05-24
	"""
	
	# libraries & modules
	import numpy as np, quaternion
	from microTools import getSCquats
	from astropy.time import Time
	from astropy.coordinates import get_body
	
	# make quaternion array from SC latitude and longitude
	lon_sc_rad = data['lon']*np.pi/180
	lat_sc_rad = data['lat']*np.pi/180
	n = np.vstack((np.zeros(np.shape(lat_sc_rad)),np.cos(lat_sc_rad)*np.cos(lon_sc_rad),np.cos(lat_sc_rad)*np.sin(lon_sc_rad),np.sin(lat_sc_rad)))
	q_coord_sc = quaternion.as_quat_array(np.transpose(n))
	
	# read SC quaternion (rotate from ECI to SC)
	qr_ECI_SC = getSCquats(int(data['gps']))
	
	# perform first rotation
	q_coord_ECI = qr_ECI_SC*q_coord_sc*quaternion.np.conjugate(qr_ECI_SC)
	
	# quaternion to rotate from ECI to Sun (place +x in Sunward direction)
	# Get sun location 
	s = get_body('sun', Time(data['gps'] ,format='gps',scale='utc'))
	sun_dec_rad = s.dec.value*np.pi/180
	sun_ra_rad = s.ra.value*np.pi/180

	# unit vector in sunward direction
	usun = np.array([np.cos(sun_dec_rad)*np.cos(sun_ra_rad),np.cos(sun_dec_rad)*np.sin(sun_ra_rad),np.sin(sun_dec_rad)])
	
	# find quaternion to go between x and sunward direction
	ux = np.array([1,0,0])
	usun_x_ux = np.cross(usun,ux)
	qr_ECIx_sun = quaternion.as_quat_array([1+np.dot(ux,usun),usun_x_ux[0],usun_x_ux[1],usun_x_ux[2]])
	qr_ECIx_sun = qr_ECIx_sun/quaternion.np.abs(qr_ECIx_sun)
	
	# perform second rotation
	q_coord_sun = qr_ECIx_sun*q_coord_ECI*quaternion.np.conjugate(qr_ECIx_sun)
	
	# extract latitude and longitude in Sunward direction
	q_coord_sun_n = quaternion.as_float_array(q_coord_sun)
	lon_sun = 180/np.pi*np.arctan2(q_coord_sun_n[:,2], q_coord_sun_n[:,1])
	lat_sun = 180/np.pi*np.arctan2(q_coord_sun_n[:,3],np.sqrt(np.square(q_coord_sun_n[:,1])+np.square(q_coord_sun_n[:,2])))
	
	# add to dictionary
	data['lon_sun'] = lon_sun
	data['lat_sun'] = lat_sun
	
	# return
	return data
	
# function to make dual corner plots
def dualCorner(data1,data2,
	keys=['Ptot','lat','lon','rx','ry','rz'],
	labels = ['$P_{tot}\,[\mu N]$','$lat\,[deg]$','$lon\,[deg]$','$r_x\,[cm]$','$r_y\,[cm]$','$r_z\,[cm]$'],
	scale = [1.0e6,1.0,1.0,100.0,100.0,100.0],
	Nbins = 30):
	"""
	function to produce a 'dual corner plot': basically a corner plot for each GRS with the 
	lower corner being GRS1 and the upper corner being GRS2. Useful for comparing chains
	Arguments:
	    data1 = dictionary containing grs1 data 
	    data2 = dictionary containing grs2 data
	    key = dictionary keys corresponding to the parameters to plot
	    labels = LaTeX labels for the keys
	    scale = scale factors for the keys
	    Nbins = number of bins
	    
	Ira Thorpe
	2018-05-30
	"""
	# import
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib
	
	# get number of keys
	Nkeys = np.shape(keys)[0]

	# initialize figure
	hf=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
	kk = 0
	# loop over keys for x (rows)
	for ii in range(0,Nkeys):
		# get x data for both GRS
		x1 = data1[keys[ii]]*scale[ii]
		N1 = np.shape(x1)[0]
		x2 = data2[keys[ii]]*scale[ii]
		N2 = np.shape(x2)[0]
		# determine x bins
		xtot = np.concatenate([x1,x2])
		xbins = np.linspace(np.min(xtot),np.max(xtot),Nbins)
		xe = xbins - 0.5*(xbins[1]-xbins[0])
		xe = np.append(xe,xe[Nbins-1]+ xe[1]-xe[0])
		# loop over keys for y (columns)
		for jj in range(0,Nkeys):
			# lower corner
			if jj < ii :
				kk = kk+1
				# get y data
				y1 = data1[keys[jj]]*scale[jj]
				y2 = data2[keys[jj]]*scale[jj]
				# determine y bins
				ytot = np.concatenate([y1,y2])
				ybins = np.linspace(np.min(ytot),np.max(ytot),Nbins)
				ye = ybins - 0.5*(ybins[1]-ybins[0])
				ye = np.append(ye,ye[Nbins-1]+ ye[1]-ye[0])
				# 2D histogram and plot
				c_x1y1,x2e,y2e = np.histogram2d(x1,y1,[xbins,ybins],normed = True)
				plt.subplot(Nkeys,Nkeys,kk)
				plt.contourf(c_x1y1,extent=[y2e.min(),y2e.max(),x2e.min(),x2e.max()],cmap=matplotlib.cm.Reds)
				ax = plt.gca()
				ax.grid(color='k',linestyle='--')
			# diagonals
			elif jj == ii :
				# histograms
				c_x1x1,x1e = np.histogram(x1,xe,normed = True)
				c_x2x2,x1e = np.histogram(x2,xe,normed = True)
				# plot
				kk = kk+1
				plt.subplot(Nkeys,Nkeys,kk)
				plt.step(xbins,c_x1x1,'r')
				plt.step(xbins,c_x2x2,'b')
				ax = plt.gca()
				ax.grid(color='k',linestyle='--')
				ax.legend(['GRS1','GRS2'])
				ax.set_yticklabels([])
			# upper corner
			elif jj > ii :
				kk = kk+1
				# determine y bins
				y1 = data1[keys[jj]]*scale[jj]
				y2 = data2[keys[jj]]*scale[jj]
				ytot = np.concatenate([y1,y2])
				ybins = np.linspace(np.min(ytot),np.max(ytot),Nbins)
				ye = ybins - 0.5*(ybins[1]-ybins[0])
				ye = np.append(ye,ye[Nbins-1]+ ye[1]-ye[0])
				# 2D histogram and plot
				c_x2y2,x2e,y2e = np.histogram2d(x2,y2,[xbins,ybins],normed = True)
				plt.subplot(Nkeys,Nkeys,kk)
				plt.contourf(c_x2y2,extent=[y2e.min(),y2e.max(),x2e.min(),x2e.max()],cmap=matplotlib.cm.Blues)
				ax = plt.gca()
				ax.grid(color='k',linestyle='--')
			# assign axes labels
			if jj == 0 :
				if ii>0 :
					ax.yaxis.label.set_text(labels[ii])
				else :
					ax.set_yticklabels([])
			elif jj == Nkeys-1 :
				if ii < Nkeys-1 :
					ax.yaxis.label.set_text(labels[ii])
					ax.yaxis.set_label_position('right')
					ax.yaxis.tick_right()
				else :
					ax.set_yticklabels([])
			else :
				ax.set_yticklabels([])
			if ii == 0 :
				ax.xaxis.label.set_text(labels[jj])
				ax.xaxis.set_label_position('top')
				ax.xaxis.tick_top()
			elif ii == Nkeys-1 :
				ax.xaxis.label.set_text(labels[jj])
			else :
				ax.set_xticklabels([])
			
	return hf

def summaryString(data,keys = ['Ptot','lat','lon','rx','ry','rz'],scale = [1.0e6,1.0,1.0,100.0,100.0,100.0]):
	"""
	function to produce a string for use in a ApJ style fancy table
	"""
	import numpy as np
	import datetime
	
	p = np.zeros([np.shape(keys)[0],3])

	for idx,kk in enumerate(keys) :
		p[idx,:] = np.percentile(data[kk]*scale[idx],[50,2.75, 97.5])

	faceNames = ['+x+x','+x+y','+y+y','+y-x','-x-x','-x-y','-y-y','-y+x','+z+z','-z-z']
	cf,bf = np.histogram(data['face'],bins=np.arange(0.5,11,1),density=True)
	if np.max(cf) > 0.7 :
		faceText = faceNames[np.argmax(cf)]
	else :
		faceText = '-'
 
	if  data['skyArea'] < (0.1*41253) :
		areaText = str('{0:.0f}'.format(data['skyArea']))
		SClatText = str('{0:.0f}'.format(data['lat_c']))
		SClonText = str('{0:.0f}'.format(data['lon_c']))
		SunlatText = str('{0:.0f}'.format(data['lat_c_sun']))
		SunlonText = str('{0:.0f}'.format(data['lon_c_sun']))
	else :
		areaText = '-'
		SClatText = '-'
		SClonText = '-'
		SunlatText = '-'
		SunlonText = '-'
	
	d = datetime.datetime.fromtimestamp(data['gps']+315964783)
	printTab = {
		'date' : d.strftime('%Y-%m-%d'),
		'gps'  : data['gps'],
		'Pmed' : p[0,0],
		'PerrU': p[0,2]-p[0,0],
		'PerrL': p[0,1]-p[0,0],
		'face' : faceText,
		'area' : areaText,
		'SClat' : SClatText,
		'SClon' : SClonText,
		'Sunlat' : SunlatText,
		'Sunlon' : SunlonText}
	

	tabStr = str(('\n{0[date]:s} & ' + 
		'{0[gps]:.0f} & ' +
		'{0[Pmed]:4.1f}^{{+{0[PerrU]:.1f}}}_{{{0[PerrL]:.1f}}} & ' +
		'{0[face]:s} & ' + 
		'{0[area]:s} & ' + 
		'{0[SClat]:s} & ' +
		'{0[SClon]:s} & ' +
		'{0[Sunlat]:s} & ' +
		'{0[Sunlon]:s} \\\\').format(printTab))

	return tabStr




############# FLATTEN LPF FUNCTIONS #############
import numpy as np
import pandas as pd

### Funtion transforms around the 'origin', or the specific spacecraft coordinte, makes all y values equal, used for Flat LPF
def transform(xs, ys, xorigin, yorigin, gotovector, index):    
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
	return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

#Returns only finite values of the dataframe, gets rid of NaN from masking, Only works with newer numpy!
def getfinite(df,i):
	# Makes places where face=! i NaN
	facenumna = df.where(df['face'] == i)
	# Makes places where NaN dissapear
	facenum = facenumna.dropna(axis = 0, how = 'any')
		
	if len(facenum['face']) == 1:
		#if there is only one value in a face, just drop it. need [] or len > 1 for hist
		facenum = df[df.face != i] 	
		#drops all values where face != face , face always = face 
		index = np.asarray(list(facenum.index.values))
			
	else:
		index = np.asarray(list(facenum.index.values)) #get index 
	return facenum, index

def dictionaryToDataFrame(dictionary):
	df = pd.DataFrame(dictionary)
	# converts new naming conventions to old
	df = df.rename(columns={'rx': 'xloc', 'ry' : 'yloc', 'rz' : 'zloc'})
	return df



#### Makes the Flattened LPF, the thing that you fold up ####
# N for each side based on area given an initial N #

def flatten_LPF(dictionary, N = 50, scale = 'log', cmap = None, 
		legend = True, add_face_nums = True, colorbar = True):
	"""
	Creates a flattened version of the LPF, 2D histogram on each face indicating where 
	the impact was
	Input:
		dictionary = data instance defined in MicroTools
	
	Arguments:
		N = integer defining how many bins histogram uses, 
			the greater N is, the longer it will take to run
			default = 50

		scale = string, defines how colormap is normalized, 
			either 'log' or 'lin'
			default = 'log'

		cmap = matplotlib colormap instance,
			default is parula from colormap.py

		legend = bool of whether to create a legend that describes percentage of hits 
			on each spacecraft face
			default = True

		add_face_nums = boolean of whether to put small numbers describing the facenumber
			default = True
	"""

	import matplotlib
	import matplotlib.pyplot as plt
	import numpy as np 
	import pandas as pd
	from matplotlib import cm
	import itertools
	import matplotlib.lines as mlines

	#3D
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib
	import matplotlib.pyplot as plt

	#3D
	from mpl_toolkits.mplot3d import Axes3D
	from mpl_toolkits.mplot3d import art3d
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection

	#2D Hist
	import matplotlib.patches as patches
	from matplotlib.path import Path

	import matplotlib.colors
	from matplotlib.colors import LogNorm
	import copy

	if not cmap:
		import colormap
		cmap = colormap.parula

	#### Defines Spacecraft coordinates ###
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

	# Converts dictionary to a pandas dataframe
	df = dictionaryToDataFrame(dictionary)


	#Cycles through faces on the LPF
	faces = np.arange(0,10)
	#Facecolor of patches, currently transparent
	alpha = 0 			

	#Parameterizing Visuals
	fig = plt.figure(figsize = (10,10))              #size of figure
	ax = fig.add_subplot(1,1,1, aspect = 'equal')    #add subplot with equal axes
	ax.set_xlim(-2,2)                                #xlim
	ax.set_ylim(-2,4)                                #ylim
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	

	#Total length and width of LPF
	Ltotal = xsc[5] - xsc[0]               
	Wtotal = ysc[2] - ysc[7]                       
	
	# Bins per unit length, Bins / LPF length
	ndensity = N / Ltotal                   

	# Bins along the Y axis on top and bottom
	binsheight = int(Wtotal * ndensity)     

	#number of Bins in the y(z) direction
	binsy = int(H * N / Ltotal) 
	
	#Find Length of sides for normalizing 
	indicies = []
	for i in faces:
		# facenumna is array of ONLY hits on face i
		facenumna = df.where(df['face'] == i)

		# Makes place where face =! i dissapear
		facenum = facenumna[np.isfinite(facenumna['face'])]
		# only a list of values that are not NaN
		indicies.append(len(list(facenum.index.values)))
	
	#colormap stuff
	if scale == 'log':   
		vmax = max(indicies)
		vmin = .5
		norm = LogNorm(vmin = vmin, vmax = vmax)
		filename = 'flatLPF_log.png'

	else:
		vmin = 0
		vmax = max(indicies) / N
		norm = matplotlib.colors.Normalize(vmin, vmax)
		filename = 'flatLPF_lin.png'

	#Sets bad values (ie: log(0) to the lowest value on the map)
	my_cmap = copy.copy(matplotlib.cm.get_cmap(cmap))
	my_cmap.set_bad(my_cmap(0))

	#Facecolors
	sidecolor = '#FF8C00' #orange
	colortop = 'navy'     #navy

	count = 0

	# finds what percentage of hits are on that face
	lennums = []
	
	#Parameterizes the faces
	# Loops through Faces to create the flattened LPF
	for i in faces:
		count += 1

		facenum, index = getfinite(df,i)
		# finds what percentage of hits are on that face
		lennums.append((len(facenum['face'])) * 100.0 / len(df['face']))

		#X and Z switched
		if i == 0:               #Parameterized Correctly, Check done
			#The left most face
			z = xsc[0]
			xpush = xsc[0] - H #Places the X correctly on the flattened LPF, does not change data
			ypush = 0

			minfacex = 0 + xpush
			maxfacex = H + xpush

			minfacey = ysc[0] + ypush
			maxfacey = ysc[1] + ypush

			width = maxfacey - minfacey
	
			# Creates the bins, based on area
			binsx = int((width / Ltotal) * N)
			binsy = int(H * ndensity)
			xs = facenum['zloc'] + xpush
			ys = facenum['yloc'] + ypush

			xs = np.asarray(xs)
			ys = np.asarray(ys) 
		 
			Hist, xedges, yedges = np.histogram2d(xs,ys, bins = [binsy,binsx], 
							range = [[minfacex, maxfacex], [minfacey, maxfacey]])

			# Transform Hist about an axis
			Hist = Hist.T 

			#Makes Patch
			xyside0 = [[minfacex, minfacey], [maxfacex, minfacey], [maxfacex, maxfacey], [minfacex, maxfacey]]
			path0 = Path(xyside0)
			patch0 = patches.PathPatch(path0, facecolor=sidecolor, lw=2, alpha = alpha)
			ax.add_patch(patch0)

			#Plots Color and clips onto patch
			ax.pcolormesh(xedges, yedges, Hist, 
					norm = norm, cmap = my_cmap, 
					clip_path = patch0, clip_on = True)

		# this side is transformed like 5
		elif i == 1:
			#base vector, pointing from 1 to 2
			basevectorx = xsc[2] - xsc[1]
			basevectory = ysc[2] - ysc[1]
			basevector = [basevectorx,basevectory]
			
			#width of base
			width = np.sqrt(basevectorx ** 2 + basevectory ** 2)
	   
			xpush = 0
			ypush = 0

			minfacex = 0 
			maxfacex = width
			
			minfacey = 0 + ypush
			maxfacey = H + ypush
			 
			binsx = int(width * ndensity) 
			
			#point that plot is turning around
			xorigin = xsc[2]
			yorigin = ysc[2]
			
			#direction transforming to, unit vector 
			gotovector = [1,0]

			#data to be transformed
			xin = facenum['xloc']
			yin = facenum['yloc']

			#transform data, flattens the side so that there are no angles
			xprime, yprime = transform(xin, yin, xorigin, yorigin, gotovector, index)
			
			#transformed data, figure out why xorigin must be added
			xs = xprime - xorigin 
			ys = facenum['zloc']

			np.asarray(xs)
			np.asarray(ys) 

			#create hist and edges from transformed data
			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx,binsy],
					range = [[minfacex,maxfacex],[minfacey,maxfacey]])
			Hist = Hist.T
			
			#find angles between sides
			# vector perpendicular to the base
			perpbase = [-1 * basevector[1], basevector[0]] 
			
			vec1 = basevector
			vec2 = [-1, 0]

			# angle between sides
			theta = 1 * np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
			
			farvec_constant =  (H / (np.linalg.norm(basevector)))
			farvec = np.multiply(farvec_constant, perpbase)
			
			#creating vertecies for patch
			xmax = xsc[1]  + farvec[0]
			ymax = ysc[1]  + farvec[1]
			cornx = xsc[2] + farvec[0]
			corny = ysc[2] + farvec[1]
			xyside1 = [[xsc[1], ysc[1]], [xsc[2], ysc[2]], [cornx, corny], [xmax, ymax]]
			 
			#places patch in the right spot
			xplace = xorigin - H * np.sin(theta) 
			yplace = yorigin - H * np.cos(theta) 
			
			#patch stuff
			path1 = Path(xyside1)
			patch1 = patches.PathPatch(path1, facecolor=sidecolor, lw=2, alpha = alpha)
			ax.add_patch(patch1)
			
			#rotate hist, rotate the histogram so that the sides are angled
			x,y = DoRotation(xedges,yedges,(theta))
			ax.pcolormesh(x + xplace, y + yplace,Hist, 
					norm = norm, cmap = my_cmap, alpha = 1, 
					clip_path = patch1, clip_on = True)
			
		### Y is actually Z
		# Topmost side that is not the octogon (8)
		elif i == 2: 
			xpush = 0
			ypush = 0 
			
			minfacex = xsc[7] + xpush
			maxfacex = xsc[6] + xpush
		
			maxfacey = H + ypush
			minfacey = 0 + ypush 
			
			width = xsc[6] - xsc[7]
			
			# Bins based on area
			binsx = int((width / Ltotal) * N)

			xs = facenum['xloc'] + xpush
			ys = facenum['zloc'] + ypush

			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx, binsy], 
					range = [[minfacex,maxfacex],[minfacey,maxfacey]])
			Hist = Hist.T 
			
			#Flips histogram up and down (visual purposes)
			Hist = np.flipud(Hist)
			xyside2 = [[minfacex, ysc[2]], [maxfacex, ysc[2]], 
					[maxfacex, ysc[2] + H], [minfacex, ysc[2] + H]]
			
			# Create patch
			path2 = Path(xyside2)
			patch2 = patches.PathPatch(path2, facecolor=sidecolor, lw=2, alpha = alpha)
			ax.add_patch(patch2)

			xedges  = np.linspace(xsc[2], xsc[3], len(xedges))
			yedges = np.linspace(ysc[2], ysc[2] + H, len(yedges))
			
			#Plots the hist
			ax.pcolormesh(xedges, yedges, Hist, norm = norm, cmap = my_cmap, 
					clip_path = patch2, clip_on=True)

		#This side is rotated 
		elif i == 3:
			
			#creates the vector pointing from vertex 4 to vertex 3, the base of side 3
			basevectorx = xsc[3] - xsc[4]
			basevectory = ysc[3] - ysc[4]
			basevector = [basevectorx, basevectory]
			
			#Length of the Base
			width = np.sqrt(basevectorx **2 + basevectory **2)
			binsx =  int(width * ndensity)                 # bins based on area
		   
			#Bins are not exactly the same, but they are pretty close 
			lenbinsx = width / binsx
			lenbinsy = H / binsy

			#point that plot is turning around
			xorigin = xsc[4]          
			yorigin = ysc[4]    
			 
			maxfacex = width 
			minfacex = 0  
		
			minfacey = 0              
			maxfacey = H              
			
			#vector points towards transformation
			gotovector = [1,0]
			
			#Data to be Transformed 
			xin = facenum['xloc']
			yin = facenum['yloc'] 
		   
			#transforms data to y = yorigin 
			xprime, yprime = transform(xin, yin, xorigin, yorigin, gotovector, index)
			
			xs = xprime - xorigin 
			ys = facenum['zloc']
			
			np.asarray(xs)
			np.asarray(ys) 

			#Creates Histogram in Easy (X,Z) reference frame
			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx,binsy],
					range = [[minfacex, maxfacex],[minfacey,maxfacey]])
			Hist = Hist.T
		 
			#vector perpendicular to the base of the side 
			perpbase = [basevector[1], -1 * basevector[0]]
		   
			#Find angle between vectors 
			vec1 = basevector
			vec2 = [1, 0]
			
			#Angle between vectors, radians
			theta = 1 * np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))) 
			#print(np.degrees(theta)) 

			farvec_constant = (H / (np.linalg.norm(basevector))) 
			farvec = np.multiply(farvec_constant, perpbase)     #Unit vector point towards top corner
			
			xmax =  xsc[3] + farvec[0]                             #X position of top right
			ymax =  ysc[3] + farvec[1]                             #Y position of top right
			cornx = xsc[4] + farvec[0]                           #X position of bot right
			corny = ysc[4] + farvec[1]                           #Y position of bot right

			# Corners for patch
			xyside3 = [[cornx, corny], [xsc[4], ysc[4]], [xsc[3], ysc[3]],
						[xmax, ymax], [cornx, corny]]

			# i dont know what these numbers are but they work
			# Constants that kept appearing
			offsetx = 0.062009 
			offsety = -1 * 0.0873899

			#Trig to figure out placement on flattened LPF
			xplace = xsc[4] + H * np.sin(theta)
			yplace = ysc[4] - H * np.cos(theta)

			path3 = Path(xyside3)
			patch3 = patches.PathPatch(path3, facecolor = sidecolor, lw = 2, alpha = alpha)
			ax.add_patch(patch3)
			
			#Rotates Matrix by theta radians
			x, y = DoRotation(xedges, yedges, (-1 * theta))
			ax.pcolormesh(x + xplace,y + yplace, Hist, 
					norm = norm, cmap = my_cmap, clip_path = patch3, clip_on = True)
			
		### X is actually Z
		elif i == 4:                            # Checked, parameterized correctly
			z = xsc[5]
			xpush = xsc[5] + H
			ypush = 0 

			maxfacex = 0 + xpush
			minfacex = -1 * H + xpush
		   
			minfacey = ysc[0] + ypush
			maxfacey = ysc[1] + ypush
			
			width = maxfacey - minfacey

			# bins based on area
			binsx = int((width / Ltotal) * N)

			xs = -1 * facenum['zloc'] + xpush
			ys = facenum['yloc'] + ypush

			Hist,xedges,yedges = np.histogram2d(xs, ys, bins = [binsy, binsx], 
					range = [[minfacex, maxfacex], [minfacey, maxfacey]])
			Hist = Hist.T

			#Create patch
			xyside4 = [[minfacex, minfacey], [maxfacex, minfacey], 
						[maxfacex, maxfacey], [minfacex, maxfacey]]
			path4 = Path(xyside4)
			patch4 = patches.PathPatch(path4, facecolor=sidecolor, lw=2, alpha = alpha)
			ax.add_patch(patch4)

			ax.pcolormesh(xedges,yedges,Hist, norm = norm, cmap = my_cmap, 
					clip_path = patch4, clip_on = True)

		#This side is transformed like 1
		elif i == 5:
			
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
		
			minfacey = 0 + ypush
			maxfacey = H + ypush

			# bins based on area
			binsx = int(width * ndensity)
			 
			#point that plot is turning around
			xorigin = xsc[6]
			yorigin = ysc[6]
			
			#direction transforming to, unit vector 
			gotovector = [1, 0]

			#data to be transformed, currently dummy data
			xin = facenum['xloc']
			yin = facenum['yloc']

			#transform data
			xprime, yprime = transform(xin, yin, xorigin, yorigin, gotovector, index)
			
			#transformed data, figure out why xorigin must be added
			xs = xprime - xorigin 
			ys = facenum['zloc']

			np.asarray(xs)
			np.asarray(ys) 

			#create hist and edges from transformed data
			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx,binsy],
					range = [[minfacex, maxfacex],[minfacey, maxfacey]])
			Hist = Hist.T
			
			#find angles between sides
			# Vector perpendicular to the base
			perpbase = [-1 * basevector[1], basevector[0]]

			vec1 = basevector
			vec2 = [-1, 0]

			# Angle between sides
			theta = np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))) #angle between sides
			
			farvec_constant =  (H / (np.linalg.norm(basevector)))
			farvec = np.multiply(farvec_constant, perpbase)
			
			#creating vertecies for patch
			xmax  = xsc[6] - farvec[0]
			ymax  = ysc[6] - farvec[1]
			cornx = xsc[5] - farvec[0]
			corny = ysc[5] - farvec[1]
			xyside5 = [[xsc[6], ysc[6]], [xsc[5], ysc[5]], [cornx, corny], [xmax, ymax], [xsc[6], ysc[6]]]
			 
			#places patch in the right spot
			xplace = xorigin + H * np.sin(theta) 
			yplace = yorigin + H * np.cos(theta) 
			
			#patch stuff
			path5  = Path(xyside5)
			patch5 = patches.PathPatch(path5, facecolor=sidecolor, lw=2, alpha = alpha)
			ax.add_patch(patch5)
			
			#rotate hist
			x,y = DoRotation(-1 * xedges, -1 * yedges, (theta))
			ax.pcolormesh(x + xplace,y + yplace, Hist, 
					norm = norm,cmap = my_cmap, alpha = 1, clip_path = patch5, clip_on = True)
			
		### Y is actually Z
		elif i == 6: 
			xpush = 0
			ypush = 0 
			
			minfacex = xsc[7] + xpush
			maxfacex = xsc[6] + xpush
		
			minfacey = 0 + ypush
			maxfacey = H + ypush 
			
			width = xsc[6] - xsc[7]

			binsx = int((width / Ltotal) * N) # bins based on area
		   

			xs = facenum['xloc'] + xpush
			ys = facenum['zloc'] + ypush

			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx, binsy], 
					range = [[minfacex, maxfacex], [minfacey, maxfacey]])
			Hist = Hist.T 
			 
			xyside6 = [[minfacex, ysc[7] - H], [maxfacex, ysc[7] - H],
					[maxfacex, ysc[7]], [minfacex, ysc[7]]]

			path6 = Path(xyside6)
			patch6 = patches.PathPatch(path6, facecolor = sidecolor, lw = 2, alpha = alpha)
			
			ax.add_patch(patch6)
			
			xedges  = np.linspace(xsc[7], xsc[6], len(xedges))
			yedges = np.linspace(ysc[7] - H, ysc[7], len(yedges))
			
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
			width = np.sqrt(basevectorx **2 + basevectory **2)
			# Bins based on area
			binsx =  int(width * ndensity)
		   
			#Bins are not exactly the same, but they are pretty close 
			lenbinsx = width / binsx
			lenbinsy = H / binsy

			#point that plot is turning around
			xorigin = xsc[7]          
			yorigin = ysc[7]        
			 
			maxfacex = width 
			minfacex = 0 
		
			minfacey = 0              
			maxfacey = H              
			
			#vector points towards transformation
			gotovector = [1,0]
			
			#Data to be Transformed 
			xin = facenum['xloc']
			yin = facenum['yloc'] 
		   
			#transforms data to y = yorigin 
			xprime, yprime = transform(xin,yin,xorigin, yorigin, gotovector, index)
			
			xs = xprime - xorigin 
			ys = facenum['zloc']
			np.asarray(xs)
			np.asarray(ys) 
			
			#Creates Histogram in Easy (X,Z) reference frame
			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [binsx, binsy],
					range = [[minfacex, maxfacex], [minfacey, maxfacey]])
			Hist = Hist.T
		   
			#vector perpendicular to the base of the side 
			perpbase = [basevector[1], -1 * basevector[0]]
		   
			#Find angle between vectors 
			vec1 = basevector
			vec2 = [1, 0]
			
			#Angle between vectors, radians
			theta = 1 * np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))) 

			farvec_constant = (H / (np.linalg.norm(basevector))) 
			# Unit vector pointing towards top corner
			farvec = np.multiply(farvec_constant, perpbase)
			
			xmax =  xsc[0] - farvec[0]                             #X position of top right
			ymax =  ysc[0] - farvec[1]                             #Y position of top right
			cornx = xsc[7] - farvec[0]                           #X position of bot right
			corny = ysc[7] - farvec[1]                           #Y position of bot right
			# Corners for patch
			xyside7 = [[cornx, corny], [xsc[7], ysc[7]], 
					[xsc[0], ysc[0]], [xmax, ymax], [cornx, corny]]

			xplace = xsc[7] - H * np.sin(theta) 
			yplace = ysc[7] + H * np.cos(theta) 

			path7 = Path(xyside7)
			patch7 = patches.PathPatch(path7, facecolor = sidecolor , lw=2, alpha = alpha)
			ax.add_patch(patch7)
			
			#Rotates Matrix by theta radians
			x, y = DoRotation(xedges, -1 * yedges, (-1 * theta))
			ax.pcolormesh(x + xplace, y + yplace, Hist, 
					norm = norm, cmap = my_cmap,
					clip_path = patch7, clip_on = True)
			
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
			# Flipped y because the bottom is viewed upside down
			ys = -1 * facenum['yloc'] + ypush
			xs = np.asarray(xs)
			ys = np.asarray(ys) 

			#Create Histogram
			Hist, xedges, yedges = np.histogram2d(xs, ys, bins = [N, binsheight], 
					range = [[minfacex, maxfacex], [minfacey, maxfacey]])
			Hist = Hist.T
			
			#Creates Patch for bottom 
			xybot = []
			for i in range(len(xsc)):
				xybot.append([xsc[i] + xpush, ysc[i] + ypush])

			pathbot = Path(xybot)
			patchbot = patches.PathPatch(pathbot, facecolor=sidecolor, lw = 2, alpha = alpha)
			patchbot1 = patchbot
			ax.add_patch(patchbot) 

			#Plots the hist, and gets cropped by the octogon
			ax.pcolormesh(xedges, yedges, Hist, norm = norm, cmap = my_cmap,
					clip_path = patchbot, clip_on = True)
	

		# This is the top, keep initial conditions, x = x, y = y ... 
		elif i ==  9: 
			
			z = H   # Zposition           
			
			#To Shift graphing Position, Must shift everything
			xpush = 0	#Shift Parameter for x
			ypush = 0 	#Shift Parameter for y

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
				xytop.append([xsc[i] + xpush, ysc[i] + ypush])

			pathtop = Path(xytop)
			patchtop = patches.PathPatch(pathtop, facecolor=colortop, lw=2, alpha = alpha)
			patchtop1 = patchtop
			ax.add_patch(patchtop) 

			#Plots the hist, and gets cropped by the octogon
			plottop = ax.pcolormesh(xedges, yedges, Hist, norm = norm, cmap = my_cmap, #interpolation='nearest', origin='lower',
					clip_path = patchtop, clip_on=True)
		
			#Makes the colorbar for all graphs, normalization is the same 
			if colorbar:
				plt.colorbar(plottop) 
	
	#Labels the facenumbers, 2 and 9 are missing because they are ugly when put on 
	if add_face_nums:
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

	if legend:
		#Makes Legend for percent hit of side numbers
		textstr = 'Percent Hit \n'
		for i in range(10):
			textstr += 'Face %i = %i'%(i, lennums[i])
			if i != 9:
				textstr +='\n'
		props = dict(boxstyle='round', facecolor='grey', alpha=0.3)

		ax.text(0.05, 0.97, textstr, transform=ax.transAxes, fontsize=10,
			verticalalignment='top', bbox=props) 

	ax.set_title('LPF Micrometeroid Impact Location %s'%(scale))
	return fig
				

