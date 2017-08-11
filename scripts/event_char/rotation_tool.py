#!/usr/bin/env python3.6
#nicolepagane|rotation_tool.py

"""
This script contains functions to convert and test for rotations from the LPF SC frame into the meteor frame.
When executed, it prompts the user to either check the z-location (sun) before rotations or to go ahead and perform rotations for all the available MCMC data.

PREREQUISITE DIRECTORIES AND FILES:
        data/run_e_*[impactGPS]*
	sky_angles/allQuats.txt

"""

import os
import numpy as np
import pandas as pd
import sys 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

############################################
## read in necessary data

# read in quaternions 
print('reading in mission quaternion data')
quat = pd.read_csv('sky_angles/allQuats.txt', sep = '\s+', names = ['time', 'q1', 'q2', 'q3', 'q4'], skiprows = 3373434 )
# ** for now, only looking at second half of mission data since the first half contains
#    no previously found events and appears to be in engineering mode due to lack of period motion (ie lissajou orbit)
print('only searching second half of attitude data')
print('')


############################################
## quaternions

def quatMult(q1,q2):
   """
   quaternion multiplication
   """
   if type(q1) == np.ndarray:
      q1 = q1.tolist()
   if type(q2) == np.ndarray:
      q2 = q2.tolist()
   if type(q1[0]) == list:
      s1 = q1[3][0]
      v1 = np.asarray([q1[0][0], q1[1][0], q1[2][0]])
   else:
      s1 = q1[3]
      v1 = np.asarray([q1[0], q1[1], q1[2]])
   if type(q2[0]) == list:
      s2 = q2[3][0]
      v2 = np.asarray([q2[0][0], q2[1][0], q2[2][0]])
   else:
      s2 = q2[3]
      v2 = np.asarray([q2[0], q2[1], q2[2]])
   s = s1 * s2 - np.dot(v1, v2)
   v = s1*v2 + s2*v1 + np.cross(v1, v2)
   v = v.reshape([len(v),1])
   q1q2 = np.vstack([s,v])
   return q1q2


def quatTrans(q, vA, inv = False):
   """
   v = [q1 q2, q3, qr]
   [q1 q2, q3] = imaginary vector
   qr = real scalar

   vB = R(qA->B) * vA
   for converting from reference A to B st q = attitude quaternion and v = coordinates and R = rotation matrix
   """
   #norm = np.sqrt(np.sum(np.asarray(q)**2))
   #if np.abs(1 - norm) < 1e-3:
   #   norm = 1
   norm = 1
   q1 = q[0]/norm
   q2 = q[1]/norm
   q3 = q[2]/norm
   qr = q[3]/norm
   q = [q1,q2,q3,qr]
   if inv == True:
      q = [-q1,-q2,-q3,qr]
   else:
      q = [q1,q2,q3,qr]
   mat = Rq(q)
   vB = mat*np.asarray(vA).reshape([3,1])
   return vB


def Rq(q):
   """
   quaternion conversion matrix
   """
   mat = np.matrix(np.asarray( [
   q[3]**2  + q[0]**2 - q[1]**2 - q[2]**2, 2*q[0]*q[1] - 2*q[3]*q[2], 2*q[0]*q[2] + 2*q[3]*q[1],
   2*q[0]*q[1] + 2*q[3]*q[2], q[3]**2 - q[0]**2 + q[1]**2 - q[2]**2, 2*q[1]*q[2] - 2*q[3]*q[0],
   2*q[0]*q[2] - 2*q[3]*q[1], 2*q[1]*q[2] + 2*q[3]*q[0], q[3]**2 - q[0]**2 - q[1]**2 + q[2]**2 ] ).reshape([3, 3]))
   return mat


def findQuat(time):
   """
   find the attitude quaternion for rotation closest to the given time
   """ 
   ind = int((time - quat.time[0])/(quat.time[len(quat.time)-1]-quat.time[0])*len(quat.time))
   indlow = ind - 10
   indup = ind + 10
   while not quat.time[indlow] <= time:
      indlow = indlow - 10
   while not quat.time[indup] >= time:
      indup = indup + 10
   acc = 100
   for i in range(indlow, indup):
      val = abs( time - quat.time[i])
      if val < acc:
         acc = val
         ind = i
   return ind
 

############################################
## additional rotations

def eqToEc(v, inv = False):
   """
   write something here later if this works
   """
   eps = 23.43701*np.pi/180
   R = np.matrix(np.asarray( [
       1, 0, 0,
       0, np.cos(eps), np.sin(eps),
       0, -np.sin(eps), np.cos(eps)  ]).reshape([3,3]))
   v = np.asarray(v).reshape([len(v),1])
   if inv == True:
      vnew = R.T*v
   else:
      vnew = R*v
   return vnew


def completeRot(events, file):
   """

   """
   gps = np.median(events[file]['t'])
   coslat = np.asarray(events[file]['lat'])
   lon = np.asarray(events[file]['lon'])
   ind = findQuat(gps)
   q = [quat.q1[ind], quat.q2[ind], quat.q3[ind], quat.q4[ind]]
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
   vCart = np.asarray(vCart).T
   alllat = []
   alllon = []
   for i in range(len(vCart)):
      X,Y,Z = quatTrans(q,vCart[i])
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
      alllat.append(lat)
      alllon.append(lon)
   return alllat, alllon


############################################
## check rotations

def checkz():
   time = quat.time
   time0 = 630806413
   timeyear = 60*60*24*365.25
   timeMarch21 = 953251215 #this is 2010 at noon bc i didn't know exact best time/year
   e = 0.0167086
   a = 1 # just make the semi-major axis 1 and scale later if necessary
   alllon = []
   alllat = []
   for i in range(0,len(time), 1000):
      offset = np.pi*2*(time[i] - time0)/timeyear
      M = np.mod(offset, np.pi*2)
      E0 = M + e*np.sin(M) + (e**2)/2*np.sin(2*M) + (e**3)/6*np.sin(3*M)
      M0 = E0 - e*np.sin(E0)
      E = E0 + (M-M0)/(1-e*np.cos(E0))
      v = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
      q = [quat.q1[i], quat.q2[i], quat.q3[i], quat.q4[i]]
      X,Y,Z = quatTrans(q, [0,0,1])
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
      alllat.append(lat)
      alllon.append(lon)
   plt.hist2d(alllon, alllat, bins = 100, cmap = 'jet', norm = LogNorm())
   plt.xlim([-np.pi, np.pi])
   plt.ylim([-np.pi/2, np.pi/2])
   plt.colorbar()
   plt.title('sun true anomaly offset st t0 = March 21')
   plt.xlabel('longitude, rad')
   plt.ylabel('latitude, rad')
   plt.show()

############################################
## make executable

if __name__ == '__main__':
   # read in sky angles
   print('reading in all spacecraft sky angles')
   files = os.listdir('data')
   for i in files:
      if i.startswith('.'):
         files.remove(i)
   for i in files:
      if i.endswith('.tgz'):
         files.remove(i)
   events = {}
   for i in files:
      events[i] = {}
      for j in ['t', 'lat', 'lon']:
         events[i][j] = []
   for iter in files:
      data = pd.read_csv(''.join(['data/', iter, '/impactchain.dat']), sep = ' ', names = ['logL', 'num', 't', 'p', 'u1', 'u2', 'coslat', 'lon', 'face', 'x', 'y', 'z'])
      data = data[int(len(data)/2):len(data)]
      events[iter]['lat'].extend(data.coslat)
      events[iter]['lon'].extend(data.lon)
      time = float(iter.split('run_e_')[1])
      events[iter]['t'].extend(data.t + time)
   check = input('do you want to check the rotations before implementing for all MCMC output in data/ subdirectories? if yes, then the rotations will be performed on the LPF z-vector (sun) for the second half of the mission attitude quaternions [yes or no]: ')
   if check == 'yes':
      checkz()
      con = input('continue rotations for MCMC data? [yes or no]: ')
      if con == 'no':
         sys.exit('exiting module now')
   print('rotating SC impact locations into meteor frame')
   if not os.path.exists('sky_angles/meteor/orientation'):
      os.makedirs('sky_angles/meteor/orientation')   
   ind = 1
   for iter in files:
      lat, lon = completeRot(events, file = iter)
      with open(''.join(['sky_angles/meteor/orientation/', iter, '.dat']), 'w') as file:
         for i in range(len(lat)):
            file.write(' '.join([str(lat[i]), str(lon[i])]))
            file.write('\n')
      print('file', ind, 'out of', len(files), 'rotated successfully')
      ind = ind + 1 

