import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
from io import BytesIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import matplotlib.colors as colors


def error_hists(pathto, currun,savedirrun,savedirmom, snrpass, make_peak, make_amp, make_q):
	
	# Initialized Figures for plotting	
	if make_peak:
		fig = plt.figure(figsize = (6,12))#, axes = plt.subplots(6,2)
		fig.suptitle('Peak Noise, Time = %s'%(currun))
	if make_q:
		figq = plt.figure(figsize = (6,12))
		figq.suptitle('Frequency vs Width, Time = %s'%(currun))
	if make_amp:
		figa = plt.figure(figsize = (6,12))
		figa.suptitle('Frequency vs Amplitude, Time = %s'%(currun))
	
	#Some constants for plotting
	d = 0
	qindex = 1  
	
	#Initializing Arrays
	params  = []
	weights = []
	
	#Where the figures are saved
	savedirrun += '/hists'
	savedirmom += '/hists'
	### If the directories don't exist, make them
	if not os.path.exists(savedirrun):
               os.mkdir(savedirrun)
	if not os.path.exists(savedirmom):
		os.mkdir(savedirmom)

	#Degrees of Freedom
	channels = [0,1,2,3,4,5]
	name_channel = ['Snx','Sny', 'Snz', 'Sntheta', 'Sneta','Snphi']
	
	#Goes through channels one at a time
	for c in channels:

		#Where the data is stored 
		freqname = pathto + '/run_e_%s/linechain_channel%s.dat'%(currun,c)
		print(name_channel[c])

		#Initializing Arrays
		freq = []  # Frequency
		amp  = []  # Amplitude
		q    = []  # Width
		
		#Opens file and goes through Line by Line
		with open(freqname) as f:
			for line in f:
				n = int(line[0])

				#If the line begins with 0, skip the line
				if n ==0:
					pass
				else:
					#Splt the line after every space
					liststrings = line.split() 
					t = 0
					
					#Adds each 3rd to the list; files are set up as n F1,A1,Q1: F2,A2,Q2		
					for i in np.arange(0,n):
						freq.append(float(liststrings[t+1]))
						amp.append(float(liststrings[t+2]))
						q.append(float(liststrings[t+3]))	
						q_index = t+3
						t = q_index
		
		#Binnumber is sqrt(length of Freq) or 1
		binnum = max(int(np.sqrt(len(freq))),1)
		bins = np.linspace(0,1,int(binnum))
		binwidth = bins[1] - bins[0]

		#Save Original Freq, Amp, Q
		freq = np.asarray(freq) 	
		amp = np.asarray(amp)	
		q = np.asarray(q)

		freqorig = freq
		amporig = amp
		qorig = q
	
		#Length Used for normalizing
		length = len(freq)
		makehist = True
		i = 0 

		#Initialize Arrays
		freqout = []     #Becomes the signal frequency
		maxfreq = []     #List of Maximum Frequencies, floats
		qout = []	 #Width of signal Frequencies
		ampout = []	 #Amplitude of signal Frequencies
		snrs  = []       #List of Signal to Noise ratios
		list_peak_freqs = [] #Full list of signal frequencies,arrays, peak1, peak2 ...
		base = 1	 #Length of base, starts as 1

		while makehist:
			#The first time through, don't mask the frequency
			if i == 0:
				freq = freq

			#If a signal is found, mask frequency (and friends)
			else:		
				#Flatten places where there is a signal
				flatsignal = np.ma.ravel(freq[signal]) 
				qsignal = np.ma.ravel(q[signal])       
				ampsignal = np.ma.ravel(amp[signal])
				
				#add flattened array to signal output
				freqout.extend(flatsignal)	
				ampout.extend(ampsignal)
				qout.extend(qsignal)

				#Append other Lists
				list_peak_freqs.append(freq[signal])
				maxfreq.append(bins[wheremax])
				
				#Redefine Freq and Friends and find next signal
				freq = freq[noise]	
				q = q[noise]
				amp = amp[noise]
			
			#Histogram the Frequency	
			n,bins = np.histogram(freq, bins = bins)	
			
			#Find the index of the largest bin	
			wheremax = np.argmax(n)	
			
			#Failed attempt at 90% confidence
			nwidth = np.where(n < n[wheremax]*.9)
			
			#width of interval to grab
			signalwidth = binwidth 
			#print('signalwidth = %s'%(signalwidth))	
		
			#Find freqs inside and outside the signal	
			noise = np.where(np.abs(freq-bins[wheremax])  > signalwidth) #outside
			signal = np.where(np.abs(freq-bins[wheremax]) < signalwidth) #inside

			#Takes the width of signal away from the base calculation, used for SNR calc
			base += -signalwidth
		
			#SNR is like a density calculation. Density of signal/density of noise	
			snr = (((len(freq[signal])/signalwidth))/(len(freq[noise])/base))			
			snrs.append(snr)
			#print('snr = %s'%(snr))

			#print(i) #Shows how many times a signal was grabbed	
			i += 1	
			
			#If SNR fails to pass, stop making the hist
			if snr < snrpass:	
				makehist = False
		j = 1
		
		#If the array is empty, make the arrays empty
		if len(freqout) < 2:
			#print('none')
			outweights = [] 
		else:
			outweights = np.ones_like(freqout)/float(length)
	
		#Weights the left over signals, (Shouldn't be empty)	
		weightnoise = np.ones_like(freq)/float(length)

		#Makes Plots of the grabbed signals and the ones left over
		if make_peak:
			ax = fig.add_subplot(6,2,d+1)
			ax.hist(freqout,bins = bins, color = '#D55B3E', weights = outweights, alpha = 1)
			ax.set_xlabel('Freq Hz')
			ax.set_ylabel('Norm Counts')
			ax.set_title('%s Peak Noise'%(name_channel[c]),size = 'small')	
			ax.set_yscale('log')
			
			ax = fig.add_subplot(6,2,d+2)	
			ax.hist(freq,bins = bins, color = '#D55B3E', weights = weightnoise, alpha = 1)
			ax.set_xlabel('Freq Hz')	
			ax.set_ylabel('Norm Counts')
			ax.set_yscale('log')
			ax.set_title('%s Noise - Peak \n SNR = %s'%(name_channel[c],snr),size='small')
		
		#Unravel arrays	
		freqout = np.ravel(freqout)
		outweights = np.ravel(outweights)
		params.append(freqout)	   #Appends all freq signals together
		weights.append(outweights) #Appends all weights togther
		

		# Makes Frequency vs Width
		if make_q:
			axq = figq.add_subplot(6,2,d+1)	
			axq.hexbin(freqorig, qorig, bins = 'log', yscale = 'log',extent = (0,1,2,8), cmap ='inferno' )
			axq.set_xlabel('Freq')	
			axq.set_ylabel('Width')
			axq.set_title('%s Total Freq vs Width'%(name_channel[c]), size='small')
			axq = figq.add_subplot(6,2,d+2)
			hist2d = axq.hexbin(freqout, qout, bins = 'log', yscale = 'log',extent = (0,1,2,8), cmap ='inferno' )
			figq.colorbar(hist2d,ax = axq)
			axq.set_xlabel('Freq')	
			axq.set_ylabel('Width')
			axq.set_title('%s Width of Peak'%(name_channel[c]), size='small')
		
		#Makes Frequency vs Amplitude	
		if make_amp:	
			axa = figa.add_subplot(6,2,d+1)	
			axa.hexbin(freqorig, amporig, bins = 'log', yscale = 'log',extent = (0,1,-21,0), cmap ='inferno' ) 
			axa.set_xlabel('Freq')	
			axa.set_ylabel('Amplitude')
			axa.set_title('%s Total Freq vs Amplitude'%(name_channel[c]), size='small')
			axa = figa.add_subplot(6,2,d+2)	
			hist2d = axa.hexbin(freqout, ampout, bins = 'log', yscale = 'log',extent = (0,1,-21,0), cmap ='inferno' )
			figa.colorbar(hist2d,ax = axa)
			axa.set_xlabel('Freq')	
			axa.set_ylabel('Amplitude')
			axa.set_title('%s Amplitude of Peak'%(name_channel[c]), size='small')
		d += 2
	
	#Saves Figures
	if make_peak:
		#Pushes plot down in order to fit suptitle
		fig.tight_layout(rect=[0,0.03,1,0.95])
		fig.savefig(savedirrun+'/params_hist.png')#%(name_channel[c]))#%(i-1)
		#fig.savefig(savedirmom+'/params_hist.png')#%(name_channel[c]))#%(i-1))
	if make_q:	
		#Pushes plot down in order to fit suptitle
		figq.tight_layout(rect=[0,0.03,1,0.95])
		figq.savefig(savedirmom+'/q_hist.png')#%(name_channel[c]))#%(i-1))
		figq.savefig(savedirrun+'/q_hist.png')#%(name_channel[c]))#%(i-1))
	if make_amp:
		#Pushes plot down in order to fit suptitle
		figa.tight_layout(rect=[0,0.03,1,0.95])
		figa.savefig(savedirmom+'/amp_hist.png')#%(name_channel[c]))#%(i-1))
		figa.savefig(savedirrun+'/amp_hist.png')#%(name_channel[c]))#%(i-1))
	

	plt.close('all')
	return params, weights

#name_channel = ['Snx','Sny', 'Snz', 'Sntheta', 'Sneta','Snphi']
#       
#       
#pathto = '/home/sophie.hourihane/website/runs'
#
##currun = '1149474616'
#
#currun = '1161032941'
#freqname = '/home/sophie.hourihane/website/runs/run_e_%s/linechain_channel5.dat'%(currun)
##freqname = '/home/sophie.hourihane/website/runs/run_e_%s/linechain_channel5.dat'%(currun)
#savedirrun = '/home/sophie.hourihane/public_html/website/public_html/images/runs/run_e_%s'%(currun)
#savedirmom = '/home/sophie.hourihane/public_html/website/public_html/images/momenta/mom_1.143809e-06'
#
#error_hists(pathto, currun, savedirrun,savedirmom, snrpass = 15)
