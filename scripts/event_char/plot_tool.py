#!/usr/bin/env python3.6
#nicolepagane|plot_tool.py 

"""
This file contains functions to generate skymaps of 2D histograms of the impact sky locations and some other random plotting tools
WRITE FUNCTION TO LOOK AT CORNER PLOT OF SKY ANGLES AND MOMENTUM AND RATE BUT FIRST FIX INTERPOLATION SPEED

PREREQUISITE DIRECTORIES AND FILES:
	data/run_e_*[impactGPS]*
	sky_angles/meteor/orientation/*[impactGPS]*
	*models/rates/[JFC or HTC]/*[impactGPS]*     *not yet bc i cant get this to run quicker

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap
import os
import sys
from catalog_tool import CI

############################################
## read in necessary data

print('reading in SC and rotated sky angles')
files = os.listdir('sky_angles/meteor/orientation/')
for i in files:
   if not i.endswith('.dat'):
      files.remove(i)

angles = {}
for i in files:
   angles[i] = {}

for i in files:
   for j in ['lat', 'lon', 'rotlat', 'rotlon']:
      angles[i][j] = []

length = []
for iter in files:
   data = pd.read_csv(''.join(['sky_angles/meteor/orientation/', iter]), sep = ' ', names = ['lat', 'lon'])
   angles[iter]['rotlat'].extend(data.lat)
   angles[iter]['rotlon'].extend(data.lon)
   data = pd.read_csv(''.join(['data/', iter.split('.dat')[0], '/impactchain.dat']), sep = ' ', names = ['logL', 'num', 't', 'p', 'u1', 'u2', 'lat', 'lon', 'face', 'x', 'y', 'z'])
   data = data[int(len(data)/2):len(data)]
   length.append(len(data))
   angles[iter]['lat'].extend(data.lat)
   angles[iter]['lon'].extend(data.lon)

print('reading in rate data for both models')


############################################
## make skymaps

def skymap(angles, indiv = False, series = False, weighted = 'length'):
   """
   this function plots the 2d histogram of the sky locations (in both the SC and meteor frame) on a mollweide projection 
   the default of the function will only plot all the events normalized by the lengths of the Markov chains
   if indiv = True, each event will be plotted individually 
   if series = True, the MCMC data will be successively appended by subsequent Markov chains to show the cumulative heatmap 
	the MCMCs are appended in order of increasing ellipsoid sky error area 
   angles is read in with the module and contains the SC and meteor event locations
   weighted = 'length' or 'ellipse'
	weights the chains by either the inverse of their length (to ensure that they are normalized) or by the inverse of their ellipsoid area (along with the normalization by length) ~the latter emphasizes better localized events
   """
   files = []
   files.extend(angles.keys())
   if not os.path.exists('sky_angles/meteor/figs'):
      os.makedirs('sky_angles/meteor/figs')
   if not os.path.exists('sky_angles/SC/figs'):
      os.makedirs('sky_angles/SC/figs')
   if indiv == True:
      if not os.path.exists('sky_angles/meteor/figs'):
         os.makedirs('sky_angles/meteor/figs')
      if not os.path.exists('sky_angles/SC/figs'):
         os.makedirs('sky_angles/SC/figs')
   if series == True:
      if not os.path.exists('sky_angles/meteor/figs/series'):
         os.makedirs('sky_angles/meteor/figs/series')
      if not os.path.exists('sky_angles/SC/figs/series'):
         os.makedirs('sky_angles/SC/figs/series')     
   low = {}
   mid = {}
   up = {}
   for i in ['lat', 'lon', 'rotlat', 'rotlon']:
      low[i] = []
      mid[i] = []
      up[i] = []
   for k in files:
      l,m,u = CI(angles[k], conf = .5)
      low['lat'].append(l[0])
      low['lon'].append(l[1])
      low['rotlat'].append(l[2])
      low['rotlon'].append(l[3])
      up['lat'].append(u[0])
      up['lon'].append(u[1])
      up['rotlat'].append(u[2])
      up['rotlon'].append(u[3])
      mid['lat'].append(m[0])
      mid['lon'].append(m[1])
      mid['rotlat'].append(m[2])
      mid['rotlon'].append(m[3])
   errar = []
   for i in range(len(files)):
      laterr = up['lat'][i] - low['lat'][i]
      lonerr = up['lon'][i] - low['lon'][i]
      errar.append(np.pi*lonerr/2*laterr/2)
   keys = files
   lat = []
   lon = []
   rotlat = []
   rotlon = []
   lbs = []
   plt.close('all')
   latbins = np.linspace(-90, 90, 100 + 1)
   lonbins = np.linspace(-180, 180, 100 + 1)
   pers = np.linspace(0, len(files)-1, len(files))
   mind = 1
   scind = 1
   for per in pers:
      errarsort = sorted(errar)
      errarsort = errarsort[int(per)]
      ind = errar.index(errarsort)
      if weighted == 'ellipse':
         lbs.extend(np.ones(length[ind])/(length[ind]*errarsort))
      else:
         lbs.extend(np.ones(length[ind])/length[ind])
      time = keys[ind].split('.dat')[0]
      lat.extend(angles[keys[ind]]['lat'])
      lon.extend(angles[keys[ind]]['lon'])
      rotlat.extend(angles[keys[ind]]['rotlat'])
      rotlon.extend(angles[keys[ind]]['rotlon'])
      if indiv == True:
         # spacecraft
         latwrap = np.arccos(np.asarray(angles[keys[ind]]['lat']))*180/np.pi
         latwrap = latwrap - 180  
         lonwrap = np.asarray(angles[keys[ind]]['lon'])*180/np.pi
         for i in range(len(lonwrap)):
            if lonwrap[i] > 180:
               lonwrap[i] = lonwrap[i] - 360
         density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True)
         plt.figure(figsize=(10,6))
         m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
         m.drawmapboundary(fill_color = 'w')
         m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
         meridian = np.arange(-180,210,30)
         m.drawmeridians(meridian, color = 'k')
         for i in range(len(meridian)):
            plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
         lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
         xs, ys = m(lonbins2d, latbins2d)
         plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
         plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
         plt.title(time)
         plt.savefig(''.join(['sky_angles/SC/figs/', time, '.jpg']))
         plt.close()
         # meteor (plz ignore how ugly all this is and how the loops [or lack thereof] are hideous)
         latwrap = np.asarray(angles[keys[ind]]['rotlat'])*180/np.pi
         for i in range(len(latwrap)):
            if latwrap[i] > 90:
               latwrap[i] = latwrap[i] - 180
         lonwrap = np.asarray(angles[keys[ind]]['rotlon'])*180/np.pi
         for i in range(len(lonwrap)):
            if lonwrap[i] > 180:
               lonwrap[i] = lonwrap[i] - 360
         plt.figure(figsize=(10,6))
         density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True)
         m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
         m.drawmapboundary(fill_color = 'w')
         m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
         meridian = np.arange(-180,210,30)
         m.drawmeridians(meridian, color = 'k')
         for i in range(len(meridian)):
            plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
         lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
         xs, ys = m(lonbins2d, latbins2d)
         plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
         plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
         plt.title(time)
         plt.savefig(''.join(['sky_angles/meteor/figs/', time, '.jpg']))
         plt.close()
      if series == True: 
         # spacecraft 
         latwrap = np.arccos(np.asarray(lat))*180/np.pi
         latwrap = latwrap - 180
         lonwrap = np.asarray(lon)*180/np.pi
         for i in range(len(lonwrap)):
            if lonwrap[i] > 180:
               lonwrap[i] = lonwrap[i] - 360
         density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True, weights = lbs)
         plt.figure(figsize=(10,6))
         m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
         m.drawmapboundary(fill_color = 'w')
         m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
         meridian = np.arange(-180,210,30)
         m.drawmeridians(meridian, color = 'k')
         for i in range(len(meridian)):
            plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
         lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
         xs, ys = m(lonbins2d, latbins2d)
         plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
         plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
         plt.title(''.join(['current ellipsoid error: ', str(errarsort)]))
         plt.savefig(''.join(['sky_angles/SC/figs/series/', str(scind), '_', time, '.jpg']))
         scind = scind + 1
         plt.close()
         # meteor (this is honestly one of the ugliest functions ive ever written, my apologies)
         latwrap = np.asarray(rotlat)*180/np.pi
         for i in range(len(latwrap)):
            if latwrap[i] > 90:
               latwrap[i] = latwrap[i] - 180
         lonwrap = np.asarray(rotlon)*180/np.pi
         for i in range(len(lonwrap)):
            if lonwrap[i] > 180:
               lonwrap[i] = lonwrap[i] - 360
         plt.figure(figsize=(10,6))
         density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True, weights = lbs)
         m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
         m.drawmapboundary(fill_color = 'w')
         m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
         meridian = np.arange(-180,210,30)
         m.drawmeridians(meridian, color = 'k')
         for i in range(len(meridian)):
            plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
         lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
         xs, ys = m(lonbins2d, latbins2d)
         plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
         plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
         plt.title(''.join(['current ellipsoid error: ', str(errarsort)]))
         plt.savefig(''.join(['sky_angles/meteor/figs/series/', str(mind), '_', time, '.jpg']))
         mind = mind + 1
         plt.close()
   # spacecraft 
   lat = np.arccos(np.asarray(lat))
   lon = np.asarray(lon)
   latwrap = lat*180/np.pi
   lonwrap = lon*180/np.pi
   for i in range(len(lonwrap)):
      if lonwrap[i] > 180:
         lonwrap[i] = lonwrap[i] - 360
   for i in range(len(latwrap)):
      if latwrap[i] > 90:
         latwrap[i] = latwrap[i] - 180
   density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True, weights = lbs)
   plt.figure(figsize=(10,6))
   m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
   m.drawmapboundary(fill_color = 'w')
   m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
   meridian = np.arange(-180,210,30)
   m.drawmeridians(meridian, color = 'k')
   for i in range(len(meridian)):
      plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
   lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
   xs, ys = m(lonbins2d, latbins2d)
   plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
   plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
   plt.title('impact sky localizations in SC frame')
   plt.savefig('SCskymap.jpg')
   plt.close()
   # meteor frame (again so sorry how ugly this is)
   lat = np.asarray(rotlat)
   lon = np.asarray(rotlon)
   latwrap = lat*180/np.pi
   lonwrap = lon*180/np.pi
   for i in range(len(lonwrap)):
      if lonwrap[i] > 180:
         lonwrap[i] = lonwrap[i] - 360
   for i in range(len(latwrap)):
      if latwrap[i] > 90:
         latwrap[i] = latwrap[i] - 180
   density, null, null = np.histogram2d(latwrap, lonwrap, [latbins, lonbins], normed = True, weights = lbs)
   plt.figure(figsize=(10,6))
   m = Basemap(projection='moll', lon_0=0, lat_0=0, resolution = 'c')
   m.drawmapboundary(fill_color = 'w')
   m.drawparallels(np.arange(-90.,120.,30.), labels = [1,0,0,0], color = 'k')
   meridian = np.arange(-180,210,30)
   m.drawmeridians(meridian, color = 'k')
   for i in range(len(meridian)):
      plt.annotate(np.str(meridian[i]),xy=m(meridian[i],-60), color = 'k')
   lonbins2d, latbins2d = np.meshgrid(lonbins, latbins)
   xs, ys = m(lonbins2d, latbins2d)
   plt.pcolormesh(xs,ys, density, cmap = 'jet', norm = LogNorm())
   plt.colorbar(orientation='horizontal', shrink=0.625, pad=0.02)
   plt.title('impact sky localizations in SC frame')
   plt.savefig('meteorskymap.jpg')
   plt.close()


############################################
## MAKE CORNER PLOTS HERE AFTER FIXING INTERPOLATION TOOL


############################################
## random functions that may or may not be useful

def impactprop(filename):
   """
   find the impact proportion from the data
   ***
   all the proprotions of impacts are expectedly 1 but this function is just to ensure so
   """
   with open('/'.join(['data', filename, 'impactchain.dat']), 'r') as file:
      data = file.readlines()
   #cut data size in half to allow for burn in
   data = data[int(len(data)/2):len(data)]
   imp = []
   for j in range(len(data)):
      data[j] = data[j].strip()
      data[j] = data[j].split()
      imp.append(int(data[j][1]))
   prop = np.sum(imp)/len(imp)
   return prop

def BF(filename):
   """
   plot bayesian factors for each event
   """
   with open('/'.join(['data', filename, 'logLchain.dat']), 'r') as file:
      data = file.readlines()
   bayes = []
   for i in range(len(data)):
      data[i] = data[i].strip()
      data[i] = data[i].split()
      bayes.append(data[i][0])
   bayes.sort()
   plt.plot(bayes)


############################################
## make executable

if __name__ == '__main__':
   indiv = input('do you want to make 2D histograms of the sky locations for every event in the SC frame and meteor frame? [yes or no]: ')
   if indiv == 'yes':
      indiv = True
   else:
      indiv = False
   print('')
   series = input('do you want to make successive histograms for all the data in order of increasing uncertainty in localization? [yes or no]: ')
   if series == 'yes':
      series = True
   else:
      series = False
   weighted = input('do you want to normalize by just the length of each chain or also by the inverse of the ellipsoid error on location such that better localized events are more prominent? [length or ellipse]: ')
   if not weighted == 'ellipse':
      weighted = 'length'
   skymap(angles, indiv = indiv, series = series, weighted = weighted)
