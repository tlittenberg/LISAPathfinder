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
    pickle.dump(data,open(str(os.cwd)+'/' + outDir+'/'+str(int(gpsTime))+'_grs' + str(int(grs)) + '.pickle','wb'))

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