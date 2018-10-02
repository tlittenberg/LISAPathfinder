# Script to read raw chain files and export them as python dictionaries
# Really only useful on Tsankawi
# Ira Thorpe
# 2018-05-18

# import libraries
import os
import fnmatch
#from pathlib import PurePath
import tarfile
from microTools import readRawChain
import shutil
import datetime

import numpy as np
import re
import pandas as pd
import numpy as np

import os
import pickle
import string



def getGPS(regex, f):
	if re.search(regex, f):
		gps = re.search(r'\d+', f).group(0)
		return gps


## Define Directories
outDir = '/local/data/ltpda/work/micrometeoroid/pyCatalog_0726'#/pyCatalog_ALL_0720'
tmpDir = '/Home/lhea/shouriha/tmp' #'/local/data/ltpda/work/micrometeoroid/tmp'
baseName = '/local/data/ltpda/work/micrometeoroid/Aug2017results'



# function to read chain files as output by MCMC tool
def readRawChain(chainDir,grs = 1, burnIn = 0.5, outDir='data'):
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

	# find directory and get gps time
	base = os.path.basename(chainDir)
	gpsTime = float(base[len(base) - 10 : len(base)])
	run = base[4]

	# load log likelihood chain
	logLfile = chainDir +'/logLchain.dat'
	dat = np.loadtxt(logLfile)
	N_full = np.shape(dat)[0]
	trim = int(float(N_full) * burnIn);
	dat = np.delete(dat, np.s_[:trim], axis=0)
	N = np.shape(dat)[0]

	# compute detection fraction N_1 / N
	N_1 = np.sum(dat[:, 0])
	dfrac =  N_1 / float(N)
	
	# load impactChain
	impFile = chainDir +'/impactchain.dat'
	dat = np.loadtxt(impFile)
	N_imp = np.shape(dat)[0]

	# We only want detections from the end, trim all but N_1
	trim = int(N_imp - N_1)
	dat = np.delete(dat, np.s_[:trim], axis=0)
	
	# build into a dictionary
	try:
		t0 = np.median(dat[:,3])
		data = {
			'segment' : gpsTime,
			'gps' : gpsTime + 1638.4 - t0,
			'N' : N,
			'dfrac': dfrac,
			'logL': dat[:, 0],
			'snr': dat[:, 1],
			't0' : -(dat[:, 3] - t0),
			'Ptot' : dat[:,4],
			'lat' : 90 - (np.arccos(dat[:,7]) * 180 / np.pi),
			'lon' : np.mod(dat[:, 8] * 180 / np.pi + 180, 360) - 180,
			'rx' : dat[:, 10],
			'ry' : dat[:, 11],
			'rz' : dat[:, 12],
			'face' : dat[:, 9],
			'run' : run
			}
	except:
		t0 = None
		data = {
			'segment' : gpsTime,
			'gps' : None,
			'N' : N,
			'dfrac': dfrac,
			'logL': None, #dat[:, 0],
			'snr':  None, #dat[:, 1],
			't0' :  None, #-(dat[:, 3] - t0),
			'Ptot' : None,# dat[:,4],
			'lat' : None, #90 - (np.arccos(dat[:,7]) * 180 / np.pi),
			'lon' : None, #np.mod(dat[:, 8] * 180 / np.pi + 180, 360) - 180,
			'rx' : None, #dat[:, 10],
			'ry' : None, #dat[:, 11],
			'rz' : None, #dat[:, 12],
			'face' : None, #dat[:, 9],
			'run' : run,  #run
			}

	# save data in processed directory
	pickle.dump(data, open(outDir + '/' + str(int(gpsTime)) + '_grs' + str(int(grs)) + '.pickle', 'wb'))
	return data


# List of Run Names
runs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 
		'h', 'i', 'j', 'k', 'l', 'm', 'n', 
		'o', 'p', 'q', 'r', 's', 't', 'u']

logFile = outDir + '/readCatalog.log'#'/local/data/ltpda/work/micrometeoroid/readCatalog.log'
log = open(logFile,'w+')
log.write('Logfile for readCatalog.py. Executed at ' + str(datetime.datetime.now())+'\n')

# loop
ii = 0
for r in runs:
	names = fnmatch.filter(os.listdir(baseName),'run_'+ r +'*.tgz')
	#log.write('Processing run ' + r + '...\n')
	for n in names:
		try:
			# Define Our regular expression to find the GPS time
			regex = r'run_' + r + r'_(\d*).tgz'
			GPS = str(getGPS(regex, n))


			# Some runs are replaced by later runs
			# we skip and do not run for those times
			isreplaced = False
			replaced = []
			if r in ['h', 'i', 'q']:
				replaced = ['r', 't', 'R', 's', 'u', 'S']
			elif r == 'k':
				replaced = ['l', 'm']
			elif r in ['a', 'b', 'c', 'd', 'e']:
				replaced = ['h', 'i']

			for rep in replaced:
				fname = baseName + '/run_' + rep + '_' + GPS + '.tgz'
				if os.path.isfile(fname):
					print('in:', r, ' replaced: ', rep)
					log.write('in: %s, replaced by: %s'%(r, rep))
					isreplaced = True
			if isreplaced:
				continue
					
			# Name of run_r_GPS file without the .tgz
			filename = 'run_' + r + '_' + GPS
			print(filename)

			# full name of .tgz files
			p = baseName + '/' + n

			# Open Tar file
			tar = tarfile.open(str(p))
			tar.extract(filename + '/' + 'impactchain.dat', path = tmpDir)
			tar.extract(filename + '/' + 'logLchain.dat', path = tmpDir)
			tar.close()

			# load impactChain
			impFile = chainDir + '/impactchain.dat'
			try:
				dat = np.loadtxt(impFile) #pd.read_csv(impFile, delimiter = '\s+')#np.loadtxt(impFile) #, delimiter = '\s')
				print(r, len(dat[0, :]))
			except IndexError or KeyError:
				print(r, 'Empty!')
			
			# Define GRS values
			if r in ['a', 'b', 'c', 'd', 'f', 'h', 
					'j', 'k', 'l', 'n', 'o', 
					'q', 'r', 't', 'R']:
				grs = 1
			elif r in ['e', 'g', 'i', 'm', 'p', 's', 'u', 'S']:
				grs = 2

			# read chain file
			readRawChain(tmpDir + '/' + filename, grs, 0.5, outDir)

			log.write('\tRead chain for ' + filename + '\n')


			# clean up file
			shutil.rmtree(tmpDir + '/' + filename)
		except:
			continue

	ii += 1

log.write('\n\nCatalog complete at ' + str(datetime.datetime.now()) + '.\n')
log.close()




