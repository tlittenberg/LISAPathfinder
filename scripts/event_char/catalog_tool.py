#!/usr/bin/env python3.6
#nicolepagane|catalog_tool.py

"""

This module contains multiple functions to read in, process, and plot micrometeorite impacts on LPF MCMC output after event reconstruction and parameter estimation.

*** NOT NECESSARY TO PLOT POSTERIOR DISTRIBUTIONS AND CREATE AN HTML VIEWER AFTER EXECUTING SOPHIE'S EVENT_VIEWER FILES ***
*** (Sophie's event_viewer scripts on github essentially do everything in this module much better and then some it where it comes to plotting and html viewers) ***
*** however this script can be run to create a catalog of impact events as a text file at 90% and 50% confidence intervals on the reconstructed parameters ***

if you would like to still make my mediocre plots and html event viewer, set the variable generateHTML to True
*** to save time and avoid redundancy, i suggest you dont do this ***

PREREQUISISTE DIRECTORIES AND FILES:
	data/run_e_*[impactGPS]*

"""

import numpy as np
import matplotlib.pyplot as plt
import corner
import math
import datetime
import os
import sys
from os.path import relpath
import re
import tarfile

############################################
## initilize catalog

def init_catalog():
   """  
   DO NOT EXECUTE THIS FUNCTION UNLESS YOU WANT TO RESTART THE CATALOG BUILDER
   if you are about to start cataloging, then execute this file to generate the catalog file

   there are only 7 parameters of interest being written into the catalog:
   impact time from the start of the segment in seconds (segment start is in the filename)
   impact momentum in N-s
   cosine latitude (in spacecraft coordinates)
   longitude (in spacecraft coordinates)
   impact location in x (m)
   impact location in y (m)
   impact location in z (m)

   """
   if os.path.exists('catalog.dat'):
      print('\nthe file already exists; you should not be executing this function.')
      ans = input('do you want to reset the catalog data (i strongly advise you not to)? [yes or no]: ')
      if ans == 'no':
         sys.exit('script terminated')
      elif ans == 'yes':
         print('the catalog file is being rewritten')
      else:
         sys.exit('script terminated. say either "yes" or "no" if this function is executed again')
   print('the catalog file is being initiated')
   with open('catalog.dat', 'w') as file:
      param = ['#gps time', 'time', 'momentum', 'cos lat [sc]', 'long [sc]', 'r', 'theta', 'phi']
      file.write('#catalog of confidence intervals (conf = .9) for parameters of interest from MCMC [low mid up]\n')
      param = '\t\t\t'.join(param)
      file.write(param)


############################################
## confidence interval

def CI(par, conf = .95):
   """
   sort the data and find the interval at the designated confidence level
   par = dictionary of parameters
   conf = confidence level
   important note: the function will not work for datasets st:
         alpha*len(data) < 1 or (conf+alpha) * len(data) > len(data) - 1
   """
   keys = []
   for i in par.keys():
      keys.append(i)
   low = []
   mid = []
   up = []
   for i in keys:
      val = sorted(par[i])
      alpha = (1 - conf)/2
      lowbound = alpha * len(val)
      upbound = (conf + alpha) * len(val)
      mid.append(np.median(val))
      if lowbound < 1 or upbound > len(val) - 1:
         print(' '.join(['Not enough data in', i]))
         low.append('NA')
         up.append('NA')
      else:
         low.append(np.asarray([val[math.ceil(lowbound - 1)], val[math.ceil(lowbound - 1)]]).mean())
         up.append(np.asarray([val[math.ceil(upbound - 1)], val[math.ceil(upbound - 1)]]).mean())
   return low, mid, up


############################################
## conversions

def convertGPS(gps):
   """
   convert GPS time to regular time. The beginning of GPS time is 6 January 1980.
   the input should be an integer/float, so ensure that the file is converted from a string
   The code in this function is inspired/modified from:
   https://stackoverflow.com/questions/33415475/how-to-get-current-date-and-time-from-event-unsegment-time-in-python
   The leap seconds were deterlowed from:
   http://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en
   gps = GPS time (as a float)
   """
   # utc = 1980-01-06UTC + (event - (leap_count(2015) - leap_count(1980)))
   utc = datetime.datetime(1980, 1, 6) + datetime.timedelta(seconds=gps - (36 - 19))
   return utc


def cartToSphere(x,y,z):
   """
   convert cartesian coordinates (x,y,z) to spherical coordinates (r, theta, phi)
   """
   x = np.asarray(x)
   y = np.asarray(y)
   z = np.asarray(z)
   r = np.sqrt(x**2 + y**2 + z**2)
   theta = np.arccos(z/np.sqrt(x**2 + y**2 + z**2))
   phi = np.arccos(x/np.sqrt(x**2 + y**2))
   return r, theta, phi

############################################
## html

def html(filename, time):
   """
   develop a html page to show the plots of the MCMC data
   filename = date/time impact of the desired plots to
   """
   with open(''.join(['./html', '/', filename, '.html']), 'w') as file:
      #make sure to include the relative path so that the webpages can be opened in a subdirectory
      img = os.listdir('/'.join([relpath('html/figs'), filename]))
      for i in img:
         if i.startswith('.') == True:
            img.remove(i)
         if i == 'figs':
            img.remove(i)
      images = []
      for i in range(len(img)):
         images.append('/'.join(['figs', filename, img[i]]))
      file.write("<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>")
      file.write('<html>')
      file.write('<head>')
      file.write(''.join(['<title>', filename, '</title>']))
      file.write('</head>')
      file.write('<body>')
      file.write(''.join(['MCMC impact data on LISA Pathfinder event: ', time]))
      file.write(''.join(['<p>', 'last edit: ', str(datetime.datetime.now().isoformat()), '</p>']))
      for i in range(len(images)):
         file.write(''.join(['<img src = ', images[i], '/>']))


def embed():
   """
   create a master html page to embed the plot pages in for easier access
   """
   #make sure to include the global path so that the webpages can be opened in a subdirectory
   link = os.listdir('/'.join([relpath('html')]))
   for i in link:
      if i.startswith('.') == True:
         link.remove(i)
      if i == 'figs':
         link.remove(i)
   hlink = []
   for i in range(len(link)):
      hlink.append('/'.join([relpath('html'), link[i]]))
   with open('catalog.dat', 'r') as file:
      data = file.readlines()
   data = data[1:len(data)]
   for i in range(len(data)):
      data[i] = data[i].strip()
      if not i == 0:
         data[i] = data[i].split()
   data[0] = data[0].strip('#')
   data[0] = data[0].split('\t')
   data[0].pop(1)
   data[0].append('')
   data = np.asarray(data).reshape([len(data), 22])
   for i in range(len(link)):
      link[i] = link[i].strip('.html')
   with open('all_events.html', 'w') as file:
      file.write("<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>")
      file.write('<html>')
      file.write('<head>')
      file.write('<title>master page</title>')
      file.write('</head>')
      file.write('<body>')
      file.write('links to plots of MCMC impact data')
      file.write('\n')
      file.write('confidence level = .90')
      file.write('\n')
      file.write('<table border = "1">')
      beige = [1,2,3,7,8,9,13,14,15,19,20,21]
      for i in range(np.shape(data)[0]):
         file.write('<tr>')
         for j in range(np.shape(data)[1]):
            if j in beige:
               file.write('<td bgcolor="LightGrey">')
            else:
               file.write('<td>')
            if i != 0 and j == 0:
               file.write(''.join(['<a href=', hlink[i-1], '>', link[i-1]]))
               file.write('</a>')
            else:
               file.write(data[i,j])
            file.write('</td>')
         file.write('</tr>')
      file.write('</body>')
      file.write('</html>')


############################################
## build catalog and make plots

def event_stat(filename, generateHTML = False):
   """
   main function to call other functions to create plots/pages and build catalog
   """
   par = {}
   #all_par = ['like', 'num_imp', 't', 'p', 'unknown1', 'unknown2', 'lat', 'lon', 'face', 'x', 'y', 'z']
   all_par = ['t', 'p', 'lat', 'lon', 'x', 'y', 'z']
   cols = [2,3,6,7,9,10,11]
   for i in all_par:
      par[i] = []
   keys = []
   for i in par.keys():
      keys.append(i)
   with open('/'.join(['data', filename, 'impactchain.dat']), 'r') as file:
      data = file.readlines()
   #cut data size in half to allow for burn in
   data = data[int(len(data)/2):len(data)]
   for j in range(len(data)):
      data[j] = data[j].strip()
      data[j] = data[j].split()
      ind = 0
      for i in keys:
         par[i].append(float(data[j][cols[ind]]))
         ind = ind + 1
   ## chanege x, y, z to r, theta, phi
   r, theta, phi = cartToSphere(par['x'], par['y'], par['z'])
   newcoord = ['r', 'theta', 'phi']
   for i in ['x', 'y', 'z']:
      par.pop(i)
   par['r'] = r
   par['theta'] = theta
   par['phi'] = phi
   ## change GPS time to regular time and make directory
   event = re.sub('run_e_', '', filename)
   dat_time = convertGPS(float(event))
   dat_time = dat_time.isoformat()
   dat_time = re.sub(':', '.', dat_time)
   print(filename, dat_time)
   ## define parameters of interest
   #plt_param = ['t', 'p', 'lat', 'lon', 'x', 'y', 'z']
   #ppar = {}
   #for i in plt_param:
   #   ppar[i] = par[i]
   #keys = []
   #for i in ppar.keys():
   #   keys.append(i)
   ppar = par
   if generateHTML == True:
      #make directory for figures
      if not os.path.exists('/'.join(['html/figs', event])):
         os.makedirs('/'.join(['html/figs', event]))
      ## plot basic statistics
      #make histograms of the parameters and show a corner plot
      gen_plot(ppar, event)
      #plot distributions and find confidence intervals (pdf, cdf)
      plot_dist(ppar, event)
   ## build up statistics by saving CI data to catalog directory
   low, mid, up = CI(ppar, conf = .90)
   with open('catalog.dat', 'a') as file:
      line = [event]
      for i in range(len(keys)):
         line.append(str('%.10e' % low[i]))
         line.append(str('%.10e' % mid[i]))
         line.append(str('%.10e' % up[i]))
      file.write('\n')
      file.write(' '.join(line))
   ## repeat building catalog but for conf = .5
   low, mid, up = CI(ppar, conf = .50)
   with open('catalog50.dat', 'a') as file:
      line = [event]
      for i in range(len(keys)):
         line.append(str('%.10e' % low[i]))
         line.append(str('%.10e' % mid[i]))
         line.append(str('%.10e' % up[i]))
      file.write(' '.join(line))
      file.write('\n')
   ## output figures on webpage
   if generateHTML == True:
      html(event, dat_time)


############################################
## plotting tools

def plot_dist(par, event):
   """
   plot the pdf and cdf for the given pareter
   par = dictionary of parameters
   event = the converted GPS date/time of the impact event
   """
   plt.close('all')
   keys = []
   for i in par.keys():
      keys.append(i)
   for i in keys:
      plt.subplot(211)
      plt.hist(par[i], bins = 100, normed = True)
      plt.title(''.join(['posterior distribution of ', i]))
      plt.subplot(212)
      plt.hist(par[i], bins = 100, cumulative = True, normed = True)
      #plt.title(''.join(['cdf of ', i]))
      plt.savefig(''.join(['html/figs/', event, '/', i, '_dist.jpg']))
      plt.clf()


def gen_plot(par, event):
   """
   make MCMC sampling and corner plots of the parameters to assess correlations and selection
   par = dictionary of parameters
   event = the converted GPS date/time of the impact event
   """
   plt.close('all')
   keys = []
   for i in par.keys():
      keys.append(i)
   iter = np.linspace(1, len(par[keys[0]]), len(par[keys[0]]))
   val = par
   for i in keys:
      plt.plot(iter, par[i])
      plt.title(''.join(['MCMC iter of ', i]))
      plt.savefig(''.join(['html/figs/', event, '/', i, 'MCMC.jpg']))
      val[i] = np.asarray(val[i]).reshape([len(val[i]), 1])
      plt.clf()
   data = val[keys[0]]
   for i in range(1, len(keys)):
      data = np.hstack([data, val[keys[i]]])
   fig = corner.corner(data, labels = keys)
   fig.savefig(''.join(['html/figs/', event, '/corner.jpg']))
   fig.clear()


############################################
## make executable

if __name__ == '__main__':
   print('')
   generateHTML = eval(input('do you want to recreate plots and html viewer (i strongly suggest no = False) [True, False]: '))
   if generateHTML == True:
      generateHTML = eval(input('are you sure? (yes = True, no = False) [True, False]: '))
   else:
      generateHTML = False
   ## check catalog file
   if not os.path.exists('catalog.dat'): 
      init_catalog()
   with open('catalog.dat', 'r') as file:
      lines = file.readlines()
   if len(lines) > 2:
      print('the catalog has', len(lines) - 2 , 'events')
   ## read in the data
   files = os.listdir('data/')
   for i in files:
      if i.startswith('.') == True:
         files.remove(i)
   tars = []
   for i in files:
      if i.endswith('.tgz') == True:
         tars.append(i)
   #untar files
   for i in tars:
      files.remove(i)
      untar = re.sub('.tgz', '', i)
      if untar not in files:
         with tarfile.open('/'.join(['data', i])) as file:
            file.extractall(path = 'data')
   files = os.listdir('./data/')
   for i in files:
      if i.startswith('.'):
         files.remove(i)
      if i.endswith('.tgz'):
         files.remove(i)
   #check to see if any of the files have been read in before and ignore if so
   if len(lines) > 2:
      new = []
      for i in range(len(files)):
         #new.append(convertGPS(float(files[i])))
         #new[i] = new[i].isoformat()
         #new[i] = re.sub(':', '.', new[i])
         new.append(files[i])
         new[i] = re.sub('run_e_', '', new[i])
      old  = []
      for i in range(2, len(lines)):
         old.append(lines[i].split()[0])
      if set(old).isdisjoint(set(new)) == False:
         incat = set(new).intersection(old)
   if 'incat' in vars():
      rep = []
      incat = tuple(incat)
      for i in incat:
         ind = new.index(i)
         new.pop(ind)
         rep.append(files.pop(ind))
      print('the file(s)', rep, 'are already in the catalog')
   if len(files) == 0:
      sys.exit('all the events in this directory are in the catalog')
   #execute tool functions to create plots and add to the catalog
   for i in files:
      event_stat(filename = i, generateHTML = generateHTML)
   if generateHTML == True:
      embed()

