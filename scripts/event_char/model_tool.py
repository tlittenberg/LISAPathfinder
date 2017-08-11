#!/usr/bin/env python3.6
#nicolepagane|model_tool.py

"""
This file contains functions to make model inferences on the micrometeorite populations (JFC = Jupiter-Family Comets, HTC = Halley-Type Comets), namely rates from extrapolated fluxes.

PREREQUISITE DIRECTORIES AND FILES:
	data/run_e_*[impactGPS]*
	sky_angles/meteor/orientation/*[impactGPS]*
	models/[JFC or HTC]_30um

"""

import numpy as np
import pandas as pd
from rotation_tool import quatTrans, eqToEc, quat
import sys
import os
import bisect 

############################################
## read in files

print('reading in JFC and HTC model data')
JFC = pd.read_csv('models/JFC_30um', sep = '\s+', names = ['lon', 'lat', 'v', 'flux'])
HTC = pd.read_csv('models/HTC_30um', sep = '\s+', names = ['lon', 'lat', 'v', 'flux'])

JFC.v = JFC.v*1000
JFC.lon = JFC.lon*np.pi/180
JFC.lat = JFC.lat*np.pi/180
HTC.v = HTC.v*1000
HTC.lon = HTC.lon*np.pi/180
HTC.lat = HTC.lat*np.pi/180
print('')

############################################
## generate extrapolated flux values

mass = {}
flux = {}
mass['JFC'] = []
mass['HTC'] = []
flux['JFC'] = []
flux['HTC'] = []

keys = []
keys.extend(flux.keys())

pgrid = np.logspace(-7, -4, 50)
key = 0
for model in [JFC, HTC]:  
   print('generating extrapolated flux values from grid of momenta values within prior range for the', keys[key], 'model')
   m = pgrid[0]/model.v*1e12
   m = np.asarray(m).reshape([len(model), 1])
   for i in range(1,len(pgrid)):
      addm = np.asarray(pgrid[i]/model.v*1e12).reshape([len(model), 1])
      m = np.hstack([m, addm])
   m = np.matrix(m)
   propm = 1800*np.power(m,-3)
   lessThan = propm > 1800*30**-3
   m[lessThan] = 0
   propm[lessThan] = 0
   m = m.T
   for i in range(len(model)):
      flux[keys[key]].append((propm[i,:]*np.matrix([model.flux[i]]*len(pgrid)).reshape([len(pgrid), 1])).tolist()[0][0])
      mass[keys[key]].append(np.median(np.asarray(m[:,i])))
   key = key + 1

print('both models with momentum dependent flux values read in')
print('')
############################################
## see new momentum v flux plots

def pVflux():
   """
   plot distribution of the model fluxes binned by momentum rather than velocity 
   """
   key = 0
   for model in [JFC, HTC]:
      pbins = np.linspace(min(mass[keys[key]])*min(model.v)*1e-12, max(mass[keys[key]])*max(model.v)*1e-12, 50)
      pmid = []
      fbins = []
      for j in range(len(pbins)-1):
         totf = 0
         for i in range(len(mass[keys[key]])):
            if pbins[j] < mass[keys[key]][i]*model.v[i]*1e-12 < pbins[j+1]:
               totf = totf + flux[keys[key]][i]
         pmid.append((pbins[j]+pbins[j+1])/2)
         fbins.append(totf)
      print(keys[key], 'fluxes binned by momentum')
      if key == 0:
         ax1 = plt.subplot(1,2,key+1)
      else:
         plt.subplot(1,2,key+1,sharex=ax1)
      plt.plot(pmid, fbins)
      plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
      plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
      plt.title(keys[key])
      key = key + 1
   plt.show()


############################################
## model selection 


def prompt_rate(mass, flux, JFC, HTC):
   """
   prompt interface to determine expected rate of impact given momentum, sky location, and time for either model
   """
   print('')
   lon = 420
   coslat = 420
   p = 420
   mod = 'blazeit'
   while not 0 <= lon <= np.pi*2:
      lon = float(input('longitude in spacecraft frame (rad) [0, 2pi]: '))
   while not -1 <= coslat <= 1:
      coslat = float(input('cosine latitude in spacecraft frame (rad) [-1, 1]: '))
   while not 1e-7 <= p <= 1e-4:
      p = float(input('momentum value (N*s) [1e-7, 1e-4]: '))
   while not mod == "JFC" and not mod == "HTC":
      mod = input('select model from which to infer rates (JFC or HTC): ')
   model = vars().get(mod)
   gps = float(input('GPS time of potential impact: '))
   if gps < quat.time[0]:
      print('')
      print('potential impact before data cutoff\n edit function such that skiprows = 31528')
   if gps > quat.time[len(quat)-1]:
      print('')
      print('ask jake to look for complete set of mission quaternions')
      print('most recent attitude data:', quat.time[len(quat)-1])
   print('')
   print('finding attitude quaternion for impact')
   ind = int((gps - quat.time[0])/(quat.time[len(quat.time)-1]-quat.time[0])*len(quat.time))
   indlow = ind - 10
   indup = ind + 10
   while not quat.time[indlow] <= gps:
      indlow = indlow - 10
   while not quat.time[indup] >= gps:
      indup = indup + 10
   acc = 100
   for i in range(indlow, indup):
      val = abs( gps - quat.time[i])
      if val < acc:
         acc = val
         ind = i
   q = [quat.q1[ind], quat.q2[ind], quat.q3[ind], quat.q4[ind]]
   # rotate spacecraft frame to meteor frame
   print('rotating sky location to meteor frame')
   time0 = 630806413
   timeyear = 60*60*24*365.25
   timeMarch21 = 953251215 #this is 2010 at noon bc i didn't know exact best time/year
   e = 0.0167086
   a = 1 # just make the semi-major axis 1 and scale later if necessary
   offset = np.pi*2*(gps - time0)/timeyear
   M = np.mod(offset, np.pi*2)
   E0 = M + e*np.sin(M) + (e**2)/2*np.sin(2*M) + (e**3)/6*np.sin(3*M)
   M0 = E0 - e*np.sin(E0)
   E = E0 + (M-M0)/(1-e*np.cos(E0))
   v = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
   vCart = [np.sin(np.arccos(coslat))*np.cos(lon), np.sin(np.arccos(coslat))*np.sin(lon), coslat]
   X,Y,Z = quatTrans(q,vCart)
   X,Y,Z = eqToEc([X.tolist()[0][0], Y.tolist()[0][0], Z.tolist()[0][0]])
   lat = Z.tolist()[0][0]
   lon = np.arctan2(Y.tolist()[0][0],X.tolist()[0][0])
   lon = lon - v + np.mod(np.pi*2*(timeMarch21-time0)/timeyear, np.pi*2)
   if lon < -np.pi:
      lon = lon + np.pi*2
   if lon > np.pi:
      lon = lon - np.pi*2
   lat = np.arccos(lat)
   lat = lat - np.pi/2
   # interpolate to infer rates
   print('interpolating momenta and sky localization to infer rate')
   #ind = int((lon - min(model.lon))/(np.pi*2)*len(model.lon))
   ind = bisect.bisect_left(model.lon, lon)
   indlow = ind - 15
   indup = ind + 15
   val = 0
   if indlow < 0:
      indlow = len(model.lon) -1 + indlow
   while not model.lon[indlow] < lon:
      indlow = indlow - 10
   while not model.lon[indup]  > lon:
      indup = indup + 10
   dist = 100
   for i in range(indlow, indup):
      latdist = model.lat[i] - lat
      londist = model.lon[i] - lon
      pdist = model.v[i]*mass[mod][i]*1e-12 - p
      distnew2 = latdist**2 + londist**2 + pdist**2
      if distnew2 < dist:
         dist = distnew2
         ind = i
   print('')
   print('nearest lon *in meteor frame* (rad) [-pi, pi]:', model.lon[ind])
   print('nearest lat *in meteor frame* (rad) [-pi/2, pi/2]:', model.lat[ind])
   print('nearest momentum value (N*s):', model.v[ind]*mass[mod][ind]*1e-12)
   print('')
   print('flux (#/m^2/s):', flux[mod][ind])
   SA = 10.64009356 # m^2
   print('rate (#/s):', flux[mod][ind]*SA)


def rate(event, time, mass, flux, JFC, HTC):
   """
   infer rates for impact events for all MCMC data 
   THIS RUNS LIKE LITERAL MOLASSES
   A BETTER SEARCH TOOL IS DEFINITELY NEEDED
   EITHER TRY BINARY SEARCH OR FIT A BETTER CURVE TO APPROX THE INITIAL INDEX GUESS
   """
   lat = event['lat']
   lon = event['lon']  
   p = event['p']
   modelind = 0
   SA = 10.64009356 # m^2
   # interpolate to infer rates
   for model in [JFC, HTC]:
      if modelind == 0: 
         mod = 'JFC'
      else:
         mod = 'HTC'
      if not os.path.exists(''.join(['models/rates/', mod])):
         os.makedirs(''.join(['models/rates/', mod]))
      with open(''.join(['models/rates/', mod, '/', time]), 'w') as file:
         for i in range(len(lon)):
            #ind = int((lon[i] - min(model.lon))/(np.pi*2)*len(model.lon))
            ind = bisect.bisect_left(model.lon, lon[i])
            indlow = ind - 15
            indup = ind + 15
            val = 0
            if indlow < 0:
               indlow = len(model.lon) -1 + indlow
               val = -np.pi*2
            while not model.lon[indlow] + val < lon[i]:
               indlow = indlow - 10
            while not model.lon[indup]  > lon[i]:
               indup = indup + 10
            dist = 100
            for j in range(indlow, indup):
               latdist = model.lat[j] - lat[i]
               londist = model.lon[j] - lon[i]
               pdist = model.v[j]*mass[mod][j]*1e-12 - p[i]
               distnew2 = latdist**2 + londist**2 + pdist**2
               if distnew2 < dist:
                  dist = distnew2
                  ind = j
            file.write(' '.join([str(flux[mod][ind]*SA), str(dist)]))
            file.write('\n')
            if np.mod(i, 50000) == 0:
               print(i/len(lon)*100, '%', 'of event', time, 'for model', mod)
      modelind = modelind + 1


############################################
## make executable

if __name__ == '__main__':
   check = input('plot the extrapolated flux values binned by momentum?\nWARNING:it takes longer than id like so maybe only do this once\nstill plot? [yes or no]: ')
   if check == 'yes':
      pVflux()
   prompt = 'yes'
   while prompt == 'yes':
      print('')
      prompt = input('do you want to infer rates for just a specific event through a text-prompted interface?\nor just determine rates for all the MCMC data? (yes = specific event, no = all MCMC) [yes or no]: ')
      if prompt == 'yes':
         prompt_rate(mass, flux, JFC, HTC)
   print('')
   check = input('infer rates for all event MCMCs? [yes or no]: ')
   if check == 'no':
      sys.exit('exiting file now')
   print('')
   # read in data
   print('reading in MCMC data')
   files = os.listdir('sky_angles/meteor/orientation/')
   for i in files:
      if not i.endswith('.dat'):
         files.remove(i)
   events = {}
   for i in files:
      events[i] = {}
      for j in ['p', 'lat', 'lon']:
         events[i][j] = []
   tind = 1
   print('')
   print('interpolating to infer rates and writing files')
   for iter in files:
      data = pd.read_csv(''.join(['sky_angles/meteor/orientation/', iter]), sep = ' ', names = ['lat', 'lon'])
      events[iter]['lat'].extend(data.lat)
      events[iter]['lon'].extend(data.lon)
      data = pd.read_csv(''.join(['data/', iter.split('.dat')[0], '/impactchain.dat']), sep = ' ', names = ['logL', 'num', 't', 'p', 'u1', 'u2', 'lat', 'lon', 'face', 'x', 'y', 'z'])
      data = data[int(len(data)/2):len(data)]
      events[iter]['p'].extend(data.p)
      rate(event = events[iter], time = iter, mass = mass, flux = flux, JFC = JFC, HTC = HTC)
      print('chain', tind, 'out of', len(files), 'completed')
      tind = tind + 1
