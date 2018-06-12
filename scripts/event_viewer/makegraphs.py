#!/usr/bin/env python

#Import the functions that I made
import lpf_graph_functions as lpf
import patch3dfunc as lpf3d
import freq_lines as freq_lines
import gif_func as gif_maker
import colormap as colormap

#Important Libraries I use

#------ Matplotlib ----#
import matplotlib
matplotlib.use('Agg')  #Backend for LIGO cluster
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.dates as mdates
import matplotlib.patches as patches			#2D Hist
from matplotlib.path import Path			#''
from mpl_toolkits.mplot3d import Axes3D			#3D Hist
from mpl_toolkits.mplot3d.art3d import Poly3DCollection #''
from mpl_toolkits.mplot3d import art3d			#''
import matplotlib.colors
from matplotlib.colors import LogNorm

#------ Other Libraries --------#
import numpy as np
import pandas as pd
import scipy as sp
import scipy.stats
import imageio

# Finding Subdirectories #
import glob
import os
import re
from os import listdir
from os.path import isfile, join

# ETC 
from decimal import Decimal
import gc
from datetime import datetime, timedelta
import copy
import argparse


#############################################Define Where the REU Folder is!!!###################################

homedir = '/home/sophie.hourihane/public_html/REU'   

####################################################################################################################

###### Mom is short for Momentum! ##############
# Linewidth 
lw = 6
#Gets subdirectories given a directory
def get_the_subdir(a_dir):
    subdir = []
    names  = []
    #print('directory = %s'%(a_dir))
    for name in os.listdir(a_dir):
        #print('name of directory = %s'%(name))
        if os.path.isdir((os.path.join(a_dir,name))):
            names.append(name)
            subdir.append((os.path.join(a_dir, name)))
    return subdir, names

#Functions transform GPS time to UTC time #
def leap(date):
    """
    Return the number of leap seconds since 6/Jan/1980
    :param date: datetime instance
    :return: leap seconds for the date (int)
    """
    if date < datetime(1981, 6, 30, 23, 59, 59):
        return 0
    leap_list = [(1981, 6, 30), (1982, 6, 30), (1983, 6, 30),
                 (1985, 6, 30), (1987, 12, 31), (1989, 12, 31),
                 (1990, 12, 31), (1992, 6, 30), (1993, 6, 30),
                 (1994, 6, 30), (1995, 12, 31), (1997, 6, 30),
                 (1998, 12, 31), (2005, 12, 31), (2008, 12, 31),
                 (2012, 6, 30), (2015, 6, 30)]
    leap_dates = map(lambda x: datetime(x[0], x[1], x[2], 23, 59, 59), leap_list)
    for j in xrange(len(leap_dates[:-1])):
        if leap_dates[j] < date < leap_dates[j + 1]:
            return j + 1
    return len(leap_dates) 
def gps2utc(week, secs):
    """
    :param week: GPS week number, i.e. 1866
    :param secs: number of seconds since the beginning of `week`
    :return: datetime instance with UTC time
    """
    secs_in_week = 604800
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    date_before_leaps = gps_epoch + timedelta(seconds=week * secs_in_week + secs)
    return date_before_leaps - timedelta(seconds=leap(date_before_leaps))

def confidence_interval(data, confidence):

    #First sorts data in order
    data = np.sort(data, axis = None)
    N = len(data)
    n_median = int(N/2)
    median = data[n_median]

    condown =((1-confidence)/2)                 #percentage down
    conup   = (confidence + ((1-confidence)/2)) #percentage up

    nup = int(N * conup)                        #index value of confidence up
    ndown = int(N * condown)                    #Index value of confidence down

    confidence_up = data[nup]
    confidence_down = data[ndown]

    return confidence_down, confidence_up, median, N

##### ----------Define Linux commands for easy use----------------#####
parser = argparse.ArgumentParser()

# Which plots to make #
parser.add_argument('-cumul','--cumulative', help = "Make the timeline (Time vs Momentum), and CDFs", action = 'store_true')
parser.add_argument('-t','--timeline', help = "Make the timeline", action = 'store_true')
parser.add_argument('--noise', help = "Make all the graphs associated with noise", action = 'store_true')
parser.add_argument('-ng','--nogifs','--nogif', help = "Make all graphs (including noise) not gifs", action = 'store_true')
parser.add_argument('-f','--fast', help = "Make the fastest plots", action = 'store_true')
parser.add_argument('-g','--gifs','--gif', help = "Make gifs", action = 'store_true')
parser.add_argument('-e','--everything', help = "Make EVERYTHING", action = 'store_true')
parser.add_argument('--old', help = "Use old Data", action = 'store_true')
parser.add_argument('-a','--skipaug', help = "Skip the runs in august", action = 'store_true')
parser.add_argument('-p','--poster', help = "Just Graphs for my Poster", action = 'store_true')
#Uses all directories by default
parser.add_argument('-d','--directory', help = "pick a run_e_directory, input gps time, (drop the run_e_)", nargs ='?', default = 'all') 
parser.add_argument('-r','--run', help = "pick a run_X directory, give the letter", nargs ='?', default = 'all') 
args = parser.parse_args()

if not args.run: #Defines where the run_e files live
	print('Using run_e')
	letter = 'e'
	where_runs = 'runs/run_e'
else:
	print 'Using run_%s'%(args.run)
	letter = args.run
	where_runs = 'runs/run_'+ args.run  #runs


#####------ Run Directories -------######
directories, names = get_the_subdir(homedir+'/' + where_runs +'/')
count = 0 #Used to print where the code is running 
count_orig = count


### Defining all Plots, Makes it easier to do switches, only have to define True #
default       = True
poster        = False			   #Makes Plots for poster

makegifs      = False			   #Makes the gif directories with images
make_all_gifs = False			   #Makes the .gif files 

timeline      = False
mom_confidence= False			   #Makes a plot of 90% confidence for all the momenta
show_path     = False

skyloc        = False			   #Makes the skylocation in spacecraft coordinates
show_path     = False			   #Makes the likelihood chain for the MCMC
mom_confidence= False			   #Makes a plot of 90% confidence for all the momenta
show_impacts  = False			   #Makes a plot of the amplitude vs time of the impact
momentum_hist = False			   #Makes momentum histogram
flatten       = False			   #Makes the Flat LPF, for making physical 3d models

confidence    = False			   #Makes a confidence interval for the noise
make_flines   = False			   #Makes hists showing where the peaks in the noise are
make_peak     = False			   #Makes a plot of where the peaks in the noise are
make_amp      = False			   #Makes a 2d hist where the peaks in amp vs freq of the noise are
make_q        = False			   #Makes a 2d hist where the peaks in the q vs freq are

print(args)
#Makes cumulative plots#

if args.poster:	
	poster        = True
	default       = False
if args.timeline:
	timeline      = True
	default       = False

if args.cumulative:			   
	if not 'all' in args.directory:
		print('### WARNING ### Cumul plots only work when all directories are selected')	
	print('Making Timeline amd Confidences')
	
	
	timeline      = True                       #Makes a timeline of all of the impacts vs momenta
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	momentum_hist = True			   #Makes momentum histogram
	default       = False
if args.fast:	
	print('Makign Fastest Plots')
	skyloc        = True			   #Makes the skylocation in spacecraft coordinates
	show_path     = True			   #Makes the likelihood chain for the MCMC
	momentum_hist = True			   ##Makes momentum histogram

	default       = False
if args.noise:
	print('Making Noise Plots')
	confidence    = True			   #Makes a confidence interval for the noise
	make_flines   = True			   #Makes hists showing where the peaks in the noise are
	make_peak     = True			   #Makes a plot of where the peaks in the noise are
	make_amp      = True			   #Makes a 2d hist where the peaks in amp vs freq of the noise are
	make_q        = True			   #Makes a 2d hist where the peaks in the q vs freq are

	default       = False
if args.nogifs:
	print('Making all Graphs, no Gifs')
	
	if not 'all' in args.directory:
		print('### WARNING ### Cumul plots only work when all directories are selected')	
	timeline      = True                       #Makes a timeline of all of the impacts vs momenta
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	
	confidence    = True			   #Makes a confidence interval for the noise
	make_flines   = True			   #Makes hists showing where the peaks in the noise are
	make_peak     = True			   #Makes a plot of where the peaks in the noise are
	make_amp      = True			   #Makes a 2d hist where the peaks in amp vs freq of the noise are
	make_q        = True			   #Makes a 2d hist where the peaks in the q vs freq are
	
	skyloc        = True			   #Makes the skylocation in spacecraft coordinates
	show_path     = True			   #Makes the likelihood chain for the MCMC
	timeline      = True
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	show_impacts  = True			   #Makes a plot of the amplitude vs time of the impact
	momentum_hist = True			   #Makes momentum histogram
	flatten       = True			   #Makes the Flat LPF, for making physical 3d models

	default       = False
if args.gifs:
	print('Making gifs')
	makegifs      = True			   #Makes the gif directories with images
	make_all_gifs = True			   #Makes the .gif files 

	default       = False


if args.everything:
	print('Making everything!')
	
	#Gifs
	makegifs      = True			   #Makes the gif directories with images
	make_all_gifs = True			   #Makes the .gif files 
	
	#Cumulative Plots
	if not 'all' in args.directory:
		print('### WARNING ### Cumul plots only work when all directories are selected')	
	timeline      = True                       #Makes a timeline of all of the impacts vs momenta
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	
	# Default plots
	skyloc        = True			   #Makes the skylocation in spacecraft coordinates
	show_path     = True			   #Makes the likelihood chain for the MCMC
	timeline      = True
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	show_impacts  = True			   #Makes a plot of the amplitude vs time of the impact
	momentum_hist = True			   #Makes momentum histogram
	flatten       = True			   #Makes the Flat LPF, for making physical 3d models
	
	#confidence    = True			   #Makes a confidence interval for the noise
	#make_flines   = True			   #Makes hists showing where the peaks in the noise are
	#make_peak     = True			   #Makes a plot of where the peaks in the noise are
	#make_amp      = True			   #Makes a 2d hist where the peaks in amp vs freq of the noise are
	#make_q        = True			   #Makes a 2d hist where the peaks in the q vs freq are

	default       = False
# Make Default Graphs, no noise or gifs 
if default:
	print('Making Default graphs')
	if not 'all' in args.directory:
		print('### WARNING ### Cumul plots only work when all directories are selected')	
	timeline      = False                      #Makes a timeline of all of the impacts vs momenta
	mom_confidence= True			   #Makes a plot of 90% confidence for all the momenta
	
	skyloc        = True			   #Makes the skylocation in spacecraft coordinates
	show_path     = True			   #Makes the likelihood chain for the MCMC
	momentum_hist = True			   #Makes momentum histogram
	show_impacts  = True			   #Makes a plot of the amplitude vs time of the impact
	flatten       = True			   #Makes the Flat LPF, for making physical 3d models
	
#Decides which directories you are using
if 'all' in args.directory:
	directories = directories
else:#lif args.directory:
	#print args.directory
	args.directory.replace('run_e_','')
	directories = [homedir+'/'+ where_runs +'/run_'+letter+'_'+ args.directory]
	#print directories
if poster:
	# Set Global fontsize 
	matplotlib.rcParams.update({'font.size': 22})
	if 'o' in args.run:
		directories = [homedir+'/' + where_runs+'/run_o_1156017899',homedir+'/' + where_runs+ '/run_o_1154963107', homedir +'/'+ where_runs +'/run_o_1157908235']
	#Old Poster
	if 'h' in args.run:
		directories = [homedir+'/' + where_runs+'/run_h_1146428350',homedir+'/' + where_runs+ '/run_h_1167944424', homedir +'/'+ where_runs +'/run_h_1144228507']
	#directories = [homedir+'/' + where_runs+'/run_e_1156209395',homedir+'/' + where_runs+ '/run_e_1154962633', homedir +'/'+ where_runs +'/run_e_1154023383']
#Initializes arrays used for Total error bar Lists

matplotlib.rcParams.update({'font.size': 22})
fSnx_list = []; fSny_list = []; fSnz_list = []; fSntheta_list=[]; fSneta_list = []; fSnphi_list = []         #Frequencies
wfSnx_list = []; wfSny_list = []; wfSnz_list = []; wfSntheta_list=[]; wfSneta_list = []; wfSnphi_list = []   #Weights (Necessary because varying lengths)
momenta_list = []
mom_weights = []
times = []
mednmomenta = []

##Initializes array for momentum confidence
mom_length   = []
mom_median   = []

mom_90_up    = []
mom_90_down  = []
mom_90_width = []

mom_50_up    = []
mom_50_down  = []
mom_50_width = []

##### Colors ########
dark_purple = '#1f2041'
light_purple = '#4b3f72'
dark_teal    = '#19647e'
light_teal   = '#119da4'
yellow = '#ffc857'
#Change how many directories you are going through if you'd like :)
#count = 22
#directories = directories[count:]
print(directories)

if timeline:
	#Gets the dates that we've looked at for the timeline
	catalogfilename = homedir + '/'+ where_runs +'/catalogfiles.txt'
	if os.path.isfile(catalogfilename):
		df_timeline = pd.read_csv(catalogfilename, header = None,delimiter = '\s+',
		names = ['times'])	
		looking_timeline = []
		
		for t in df_timeline['times']: 
			looking_time = int(t[2:12])
			len_time     = int(t[13:17])	
			looking_timeline.append([gps2utc(0,looking_time),gps2utc(0,looking_time + len_time)])

		print(looking_timeline)
for direc in directories:

	if not 'run' in direc[len(homedir+'/'+ where_runs +'/'):-13]:
		continue 
for direc in directories:
    
    #Skips directories that are not run directories
    if not 'run' in direc[len(homedir+'/'+ where_runs +'/'):-13]:
	continue 
    else: 

	if args.skipaug:	 
		rundir = direc[-16:]  #'run_e_1144228507'
		currun = rundir[-10:] #'1144228507
		print (1185580818. - float(currun))
		if float(currun) < 1156636817:
			if float(currun) > 1154044817: 
				print "This run is in August 2016"
				continue
			else:
				pass	
	print('count graph plots = %s'%(count))
	print('going through files in %s'%(direc))
	
	total_hist_path = homedir +'/'+ where_runs +'/total_hists'
	if not os.path.exists(total_hist_path):
		os.mkdir(total_hist_path)
	
	#Gets only files in the directories
	onlyfiles = [fil for fil in listdir('%s'%(direc)) if isfile(join('%s'%(direc),fil))]


	#Name of directories, will maybe have to change on different computer
	rundir = direc[-16:]  #'run_e_1144228507'
	currun = rundir[-10:] #'1144228507
	currun_title = "%s, %s"%(gps2utc(0,float(currun)),currun)	
	
	#Impact Chain Info
	impactname = direc + '/impactchain.dat'
	
	#### Import Data ####
	if 'e' in letter:
	    df = pd.read_csv(impactname, header = None,delimiter = '\s+',
	    names = ['logp','impactnum','time','mom','whatever','who??','coslat', 'longi', 'face', 'xloc','yloc','zloc'])#, skiprows = int(len(index)/2))
	else: 
	    df = pd.read_csv(impactname, header = None,delimiter = '\s+',
	    names = ['logp','SNR','impactnum','time','mom','whatever','who??','coslat', 'longi', 'face', 'xloc','yloc','zloc'])#, skiprows = int(len(index)/2))

	#Only reads second 1/2 of Data, MCMC burn in
	print("Length of dataframe = %s"%(len(df['impactnum'])))
	df = df.tail(n = int(len(df['impactnum'])/2))
	#print('DF length = %s'%(len(df['impactnum'])))

	#Find Median Momentum for title
	mednmom = np.median(df['mom'])
	mednmom = "{:.4e}".format(mednmom) #Scientific Notation, 6 digits
	print('momentum median = %s'%(mednmom))
	momdir = 'mom_%s'%(mednmom)
	rundir = 'run_e_%s'%(currun)

	#Makes Momenta and Run Directories for storing images
	
	#dirpath = '/home/sophie.hourihane/public_html/website/public_html/images/momenta/%s'%(momdir)
	dirpath = direc + '/images_%s_%s'%(currun,mednmom)
	if not os.path.exists(dirpath):
		os.mkdir(dirpath)

	rundirpath = dirpath 
	
	##### Plot confidence, Noise in Data, Miami Dolphins colors#####
	if confidence:
		lpf.ploterror(mednmom, currun, dirpath, rundirpath, homedir, where_runs) #Yscale is always Log

	####### Runs Flat LPF Functions, the ones that fold up in real life ########
	N = 50
	if flatten:
		print('Flattening LPF')
		lpf.flatten_LPF(df,N, 'lin', mednmom, currun, colormap.parula, dirpath, rundirpath)
		lpf.flatten_LPF(df,N, 'log', mednmom, currun, colormap.parula, dirpath, rundirpath)
	#Makes Timeline of the impacts
	if timeline:

		gpstime = int(float(currun))	
		utctime = gps2utc(0,gpstime)
		#print(gpstime,utctime)
		times.append(utctime)
		mednmomenta.append(float(mednmom))
	####### Momentum Hist ... Histograms of the momenta distributions ########
	if momentum_hist:
		print('Momentum Histogram') 

		fig = plt.figure(figsize = (7,6))
		ax = fig.add_subplot(111)

		purple = colormap.parula(1)#'#FFC13C'

		binnum = int(np.sqrt(len(df['mom'])))
		
		mom_wei = np.ones_like(df['mom'])/len(df['mom'])
		hist,bins,patches = ax.hist(df['mom'], bins = binnum,weights = mom_wei, color = purple)
		#print('hist sum individual = %s'%(np.sum(hist)))#/len(df['impactnum'])))
		momenta_list.extend(df['mom'])
		mom_weights.extend(mom_wei)
		
		ax.set_title('Momentum Distribution \n Median Momentum = %s kg m/s \n %s s'%(mednmom, currun_title))#, y =1.08)
		#ax.title(fig_title, y=1.08)
		ax.set_xlabel('Momentum Kg m / s')
		ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

		ax.set_ylabel('PDF')
		plt.tight_layout()#rect=[0,.03,1,0.95])
		plt.savefig('%s/momentum_%s_%s.png'%(dirpath,currun,mednmom))
		#plt.savefig('%s/momentum_%s.png'%(rundirpath,currun))
		plt.close()
		plt.close()
	####### Path: Shows how the MCMC moves through space############
	if show_path:
		print('Likelihood Chain')
		fig = plt.figure(figsize = (7,6))
		ax = fig.add_subplot(1,1,1)

		light_purple = colormap.parula(0) #'#5F2CFF'

		index = list(df['logp'].index.values)                       #Lets you skip NaN stuff
		ax.plot(index, df['logp'], color = light_purple)
		ax.set_title('Log Likelihood of Impact Chain, \n Median Momentum = %s kg m/s \n %s s'%(mednmom, currun_title), y=1.03)#, size='small')
		#ax.title(fig_title, y=1.08)
		ax.set_xlabel('Step in Chain')
		ax.set_ylabel('Log Likelihood')
		ax.set_xlim(min(index))

		plt.tight_layout()
		#plt.tight_layout()
		plt.savefig('%s/logprob_chain_%s_%s.png'%(dirpath,currun,mednmom))
		#plt.savefig('%s/logprob_chain_%s.png'%(rundirpath,currun))
		plt.close()
		plt.close()

	##### Skyloc: Plots origin of Micrometeorite in Spacecraft coordinates #########
	if skyloc:
		print('Sky Location')
		fig = plt.figure(figsize = (7,6))
		ax = fig.add_subplot(1,1,1)
		norm = LogNorm()
		my_cmap = colormap.parula
		my_cmap.set_bad(my_cmap(0))
		binnum = N	
		#Xlim
		longlim = [0,2*np.pi]
		#Ylim
		coslatlim = [-1,1]
		
		#binnumber, sqrt of length
		binshist = np.sqrt(len(df['longi']))#2000#int(2*np.pi*(binnum/(max(df['longi']) - min(df['longi']))))
		
		counts,xedges,yedges,im= ax.hist2d(df['longi'],df['coslat'], bins = binshist,#np.sqrt(len(df['longi'])),#binshist, 
				 range = np.array([longlim,coslatlim]),norm = norm,cmap = my_cmap)
		ax.set_title('Sky Position of Impact, Spacecraft Coordinates\n Median Momentum = %s \n %s'%(mednmom, currun_title))#, size ='small')
		ax.set_xlabel('Longitude')
		ax.set_ylabel('Cosine Latitude')
		ax.set_xlim(0,2*np.pi)
		ax.set_ylim(-1,1)
		plt.colorbar(im)

		plt.tight_layout()#rect=[0,.03,1,0.95])
		#plt.tight_layout()
		plt.savefig('%s/skyloc_%s_%s.png'%(dirpath,currun,mednmom))
		
		plt.close()        
		plt.close()
		plt.close('all')
	### Makes graphs of Peak Freq, Amp vs Freq, and Q vs Freq, Instrument performance over time ##
	if make_flines:
		print('Making Peak Frequency Graphs')
	
		pathto = homedir + '/' + where_runs #rundirpath.replace(rundir,'')#'/home/sophie.hourihane/website/runs'
		
		#Arbitrary Passing value for how intense the lines must be to grab. Make it > 1. 
		snrpass = 5
		
		#What directories to save images in
		savedirrun = rundirpath
		savedirmom = dirpath
		
		#Params and weight len = 6, but filled with stuff 
		par,wei = freq_lines.error_hists(pathto, currun,savedirrun,savedirmom, snrpass, make_peak, make_amp, make_q)
		#Change the lists to numpy arrays
		par = np.asarray(par)
		wei = np.asarray(wei)
		
		#Adds frequency for each parameter to a total list
		fSnx_list.extend(par[0])
		fSny_list.extend(par[1])
		fSnz_list.extend(par[2])
		fSntheta_list.extend(par[3])
		fSneta_list.extend(par[4])
		fSnphi_list.extend(par[5])	

		#Adds weights for each parameter to total lists
		wfSnx_list.extend(wei[0])
		wfSny_list.extend(wei[1])
		wfSnz_list.extend(wei[2])
		wfSntheta_list.extend(wei[3])
		wfSneta_list.extend(wei[4])
		wfSnphi_list.extend(wei[5])	
	if show_impacts:
		print('Making Impact Acceleration Model')
		#paramnum = [0,1,2,3,4,5,6]
		fig,axs = plt.subplots(nrows = 6, ncols = 1, figsize = (8,12), sharex = True)
		#plt.subplots_adjust(left=0.2, wspace=0.8, top=0.8)
		plt.suptitle('Impact Reconstruction \n Momentum = %s [Ns] \n %s'%(mednmom, currun_title +'[s]'))
		
		paramname = ['X', 'Y', 'Z', 'Theta', 'Eta', 'Phi']
		for p in range(len(paramname)):
			filename  = direc + '/TD_model_%s.dat'%(p) 
			df_model = pd.read_csv(filename, header = None,delimiter = '\s+',
			names = ['time','amp'], skiprows = 1)#, skiprows = int(len(index)/2))
			
			filename  = direc + '/TD_data_%s.dat'%(p) 
			df_data = pd.read_csv(filename, header = None,delimiter = '\s+',
			names = ['time','amp'], skiprows = 1)#, skiprows = int(len(index)/2))
			ax = axs[p]
			ax.set_ylabel(paramname[p] + ' amp')
			ax.plot(1638.4-df_data['time'], df_data['amp'], label = 'Data',color = dark_purple, zorder = 1)
			ax.plot(1638.4-df_model['time'],df_model['amp'],label ='Model', color = light_teal, zorder = 10)
			#ax.set_yscale('log')
			
		ax.set_xlabel('Time [s]')
		ax = axs[0]
		ax.legend(loc = 'upper left')
		filename = rundirpath + '/impact_pulse_%s_%s.png'%(currun,mednmom)
		plt.tight_layout(rect=[0,.03,1,0.93])
		plt.savefig(filename)
		plt.close('all')
			
	if mom_confidence:
		down90, up90, median, length = confidence_interval(df['mom'], confidence = .90)
		width_90 = up90 - down90
		mom_90_up.append(up90)
		mom_90_down.append(down90)
		mom_90_width.append(width_90)

		down50, up50, median, length = confidence_interval(df['mom'], confidence = .50)					
		width_50 = up50 - down50
		mom_50_up.append(up50)
		mom_50_down.append(down50)
		mom_50_width.append(width_50)
		
		mom_median.append(median)	
	####Make and fill image folders for gifs, making the .gif is a separate function###
	if makegifs:
		
		gifpath = dirpath + '/gifs'
		if not os.path.exists(gifpath):
			os.mkdir(gifpath)

		#arbitrary binnumber, this function is SLoW 
		N = 50

		lpf3d.patch3d_LPF(df,N,'log', mednmom, currun, gifpath, cmap = colormap.parula)
	   
	if make_all_gifs:
		gif_maker.gif_maker(gifdir = dirpath + '/gifs', currun = currun, mednmom = mednmom, savedir = dirpath)
	if poster:
		posterpath = homedir+'/poster'
		if not os.path.exists(posterpath):
			os.mkdir(posterpath)
		
		fig = plt.figure(figsize = (10,10))
		ax = fig.add_subplot(1,1,1)
		norm = LogNorm()
		my_cmap = colormap.parula
		my_cmap.set_bad(my_cmap(0))
		binnum = N	
		#Xlim
		longlim = [0,2*np.pi]
		#Ylim
		coslatlim = [-1,1]
		
		#binnumber, sqrt of length
		binshist = np.sqrt(len(df['longi']))#2000#int(2*np.pi*(binnum/(max(df['longi']) - min(df['longi']))))
		
		counts,xedges,yedges,im= ax.hist2d(df['longi'],df['coslat'], bins = binshist,#np.sqrt(len(df['longi'])),#binshist, 
				 range = np.array([longlim,coslatlim]),norm = norm,cmap = my_cmap)
		ax.set_title('Impact Sky Position')#, size ='small')
		ax.set_xlabel('Longitude')
		ax.set_ylabel('Cosine Latitude')
		ax.set_xlim(0,2*np.pi)
		ax.set_ylim(-1,1)
		plt.tight_layout()#rect=[0,.03,1,0.95])
		plt.savefig('%s/skyloc_%s_%s.png'%(posterpath,currun,mednmom))
		plt.close('all')
		
		print('Momentum Histogram') 

		fig = plt.figure(figsize = (10,10))
		ax = fig.add_subplot(111)

		purple = colormap.parula(1)#'#FFC13C'

		binnum = int(np.sqrt(len(df['mom'])))
		
		mom_wei = np.ones_like(df['mom'])/len(df['mom'])
		hist,bins,patches = ax.hist(df['mom'], bins = binnum,weights = mom_wei, color = purple)
		#print('hist sum individual = %s'%(np.sum(hist)))#/len(df['impactnum'])))
		momenta_list.extend(df['mom'])
		mom_weights.extend(mom_wei)
		
		ax.set_title('Momentum PDF')#, y =1.08)
		#ax.title(fig_title, y=1.08)
		ax.set_xlabel('P [Ns]')
		ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

		ax.set_ylabel('p(P|d)')
		plt.tight_layout()#rect=[0,.03,1,0.95])
		plt.savefig('%s/momentum_%s_%s.png'%(posterpath,currun,mednmom))
		#plt.savefig('%s/momentum_%s.png'%(rundirpath,currun))
		plt.close()
	
	del df #cleans up arrays for garbage collector
	
	#Garbage collector shows us whats's up <3 
	n = gc.collect()
	count += 1
	print('garbage = %s'%(n)) 
	

if poster:
	#paramnum = [0,1,2,3,4,5,6]
	fig,ax = plt.subplots(nrows = 1, ncols = 1, figsize = (20,5))
	#plt.subplots_adjust(left=0.2, wspace=0.8, top=0.8)
	#plt.suptitle('Impact Reconstruction'%(mednmom, currun_title +'[s]'))
		
	paramname = ['X', 'Y', 'Z', 'Theta', 'Eta', 'Phi']
	runext = '/'+ where_runs + '/run_h_1144228507'
	
	#'/run_o_1157908235' #'/run_e_1154962633'

	filename  = homedir + runext + '/TD_model_%s.dat'%(0) 
	df_model = pd.read_csv(filename, header = None,delimiter = '\s+',
	names = ['time','amp'], skiprows = 1)#, skiprows = int(len(index)/2))
	
	filename  = homedir + runext + '/TD_data_%s.dat'%(0) 
	df_data = pd.read_csv(filename, header = None,delimiter = '\s+',
	names = ['time','amp'], skiprows = 1)#, skiprows = int(len(index)/2))
	ax = ax
	ax.set_ylabel('Amp')
	ax.plot(1638.4-df_data['time'], df_data['amp'], label = 'Data',color = dark_purple, zorder = 1)
	ax.plot(1638.4-df_model['time'],df_model['amp'],label ='Model', color = light_teal, zorder = 10)
		#ax.set_yscale('log')
	ax.set_ylim(75,-75)#-30,20)#-13,13)
	ax.set_xlim(1300,1500)#1000, 1250) #750,1000)
	ax.set_xlabel('Time [s]')
	ax.set_title('Impact Reconstruction\n2016-04-09')#2016-09-14') #2016-08-11')
	
	ax.legend(loc = 'upper left')
	filename = posterpath + '/impact_pulse_%s_%s.png'%(currun,mednmom)
	plt.tight_layout()#rect=[0,.03,1,0.93])
	plt.savefig(filename)
	plt.close('all')

if mom_confidence:
	print('Making Momentum CDFs')	
	path = homedir + '/' + where_runs #/home/sophie.hourihane/public_html/website/public_html'

	########## Plots Cumulative momentum WRT median ####################
	fig = plt.figure(figsize = (10,10))
	ax = fig.add_subplot(1,1,1)
	mom_median = np.asarray(mom_median)
	mom_90_up = np.asarray(mom_90_up)
	mom_90_down = np.asarray(mom_90_down)
	mom_90_width = np.asarray(mom_90_width)
	mom_50_width = np.asarray(mom_50_width)
	mom_50_up = np.asarray(mom_50_up)
	mom_50_down = np.asarray(mom_50_down)
	med_order = mom_median.argsort() #returns ordered indicies of mom_median
	med_order = np.asarray(med_order)
	med_med   = mom_median[med_order]
	
	ranks = np.linspace(0,len(mom_median),len(mom_median))/len(mom_median)
	#print(ranks)
	ax.set_yscale('linear')
	ax.set_xscale('log')
	lw = 4
	up90 = ax.plot(mom_90_up[med_order],ranks,label = '90% Credible', color = dark_purple, linewidth = lw)	
	down90 = ax.plot(mom_90_down[med_order],ranks, color = dark_purple, linewidth = lw)
	#ax.fill_betweenx(med_med, down90,up90, color = light_purple, interpolate = True)
	
	up50 = ax.plot(mom_50_up[med_order],ranks,label = '50% Credible', color = dark_teal, linewidth = lw)	
	down50 = ax.plot(mom_50_down[med_order],ranks, color = dark_teal, linewidth = lw)	
	#ax.fill_betweenx(med_med, down50,up50, color = light_teal, interpolate = True)

	ax.plot(med_med,ranks,label = 'Median', color = yellow, linewidth = lw)
	ax.set_title('Momentum CDF')
	ax.set_xlabel('Momentum Ns')
	ax.set_ylabel('CDF')
	ax.legend(loc='upper left')
	plt.tight_layout()
	plt.savefig(path+'/total_hists/mom_confidence.png')
	plt.close()

	########### Plots Cumulative momentum WRT Fractional Error ########

	fig = plt.figure(figsize = (10,10))
	ax = fig.add_subplot(1,1,1)
	width_over_median = mom_90_width/mom_median	

	width90 = width_over_median.argsort() #returns ordered indicies of mom_median
	width90 = np.asarray(width90)
	width90_med   = mom_median[width90]
	
	ranks = np.linspace(0,len(mom_median),len(mom_median))/len(mom_median)
	#print(ranks)
	#ax.set_yscale('linear')
	#ax.set_xscale('log')
	lw = 6
	ax.plot(width_over_median[width90], ranks,label = 'Fractional Error', color =light_purple, linewidth = lw)
	ax.set_title('Fractional Error CDF')
	ax.set_xlabel('Fractional Error')
	ax.set_ylabel('CDF')
	#ax.legend(loc='upper left')
	path = homedir + '/'+ where_runs  #'/home/sophie.hourihane/public_html/website/public_html'
	plt.tight_layout()
	plt.savefig(path+'/total_hists/mom_90_confidence.png')
	plt.close()


if timeline:
	# Make the timeline 
	# Set Global fontsize 
	matplotlib.rcParams.update({'font.size': 22})
	
	fig = plt.figure(figsize = (20,3))
	frame = plt.gca()
	ax = fig.add_subplot(1,1,1)
	#ax.set_frame_on(False)
	title = ax.set_title('Impact Timeline')
	ones = np.ones(len(times))
	ax.scatter(times,ones, c=dark_purple, s=80)
	
	if os.path.isfile(catalogfilename):
		#makes a line for every second interval of time, takes a bit of time
		for t in looking_timeline:
			#print(t)
			ax.plot(t,[1.25,1.25], c = light_teal , linewidth = lw)
	ax.get_yaxis().set_ticks([])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
           #marker='s', s=100)
	

	fig.autofmt_xdate()
	ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

	path = homedir + '/'+ where_runs #'/home/sophie.hourihane/public_html/website/public_html'
	plt.tight_layout()
	plt.savefig(path+'/total_hists/timeline.png')
	plt.close()

# Makes a histogram of all the momenta. It isn't super helpful, but its here
if momentum_hist:	
	fig = plt.figure(figsize = (6,10))	
	binnum = int(np.sqrt(len(momenta_list)))
	mom_weights = np.asarray(mom_weights)/(count - count_orig)
	ax = fig.add_subplot(3,1,1)	
	title = ax.set_title('Total Momenta Hists')#, size = 'small') 
	hist,bins,patches = ax.hist(momenta_list,bins = binnum, weights = mom_weights)
	ax.set_ylabel('Counts')
	#print('hist sum= %s'%(np.sum(hist)))
	#print('hist total sum individual = %s'%(np.sum(hist)/len(momenta_list)))
	ax = fig.add_subplot(3,1,2)	
	ax.hist(momenta_list,bins = binnum, weights = mom_weights)
	ax.set_yscale('log')
	ax.set_ylabel('Counts')

	ax = fig.add_subplot(3,1,3)	
	ax.hist(momenta_list,bins = binnum, weights = mom_weights)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('Momentum kg m/s')

	ax.set_ylabel('Counts')
	path = homedir + '/' + where_runs #'/home/sophie.hourihane/public_html/website/public_html'
	plt.tight_layout()
	plt.savefig(path+'/total_hists/total_momenta.png')
	
if make_flines: #Makes hists of where the errors are for ALL TIME, probably not weighted correctly smh 
	#Unravel Frequencies in DOFs
	Snx = np.ravel(fSnx_list)
	Sny = np.ravel(fSny_list)
	Snz = np.ravel(fSnz_list)
	Sntheta = np.ravel(fSntheta_list)
	Sneta   = np.ravel(fSneta_list)
	Snphi   = np.ravel(fSnphi_list)	
	
	#Unravel Weights in DOFs
	wSnx = np.ravel(wfSnx_list)
	wSny = np.ravel(wfSny_list)
	wSnz = np.ravel(wfSnz_list)
	wSntheta = np.ravel(wfSntheta_list)
	wSneta   = np.ravel(wfSneta_list)
	wSnphi   = np.ravel(wfSnphi_list)	
	params = [Snx,Sny,Snz,Sntheta,Sneta,Snphi]
	
	weights = [wSnx,wSny,wSnz,wSntheta,wSneta,wSnphi]
	nameparams = ['Snx','Sny','Snz','Sntheta','Sneta','Snphi']
	path = homedir +'/'+ where_runs #'/home/sophie.hourihane/public_html/website/public_html'
	i = 0
	
	fig = plt.figure(figsize = (6,12))
	for p in params:
		if len(p) > 2:
			binnum = np.sqrt(len(p))
			bins = np.linspace(0,1,binnum)
		else: 
			binnum = 50
			bins = np.linspace(0,1,binnum)
		
		ax = fig.add_subplot(6,1,i+1)
		ax.hist(p,bins = bins, weights = weights[i])
		ax.set_yscale('log')
		title = ax.set_title('Freq %s'%(nameparams[i]))#, size = 'small') 
		i += 1
	plt.tight_layout()
	plt.savefig(path+'/total_hists/total_params.png')
