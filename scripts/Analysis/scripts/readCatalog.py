# import
import numpy as np, quaternion
import os
import pathlib
import pickle
import matplotlib.pyplot as plt
import matplotlib
# For printing exceptions
import traceback

import datetime
from impactClass import impactClass
from populationClass import population as pop
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-g", "--gif",
		help="default, false, decides whether or not to create gif directories", 
		action="store_true")

args = parser.parse_args()

if not args.gif:
	print('WARNING: You are not creating gifs\n \tto create: python readCatalog.py -g')


# list all of the data files and initialize the catalog LaTeX
# setup directory structure
p = pathlib.PurePath(os.getcwd())
BASE_DIR = str(p.parent)
dataDir = '/data'

### MIGHT NEED TO UPDATE THIS ###
#impactList = impactClassList(grs = 1, getValid = True, BASE_DIR = BASE_DIR, directory = '/data/ONLY_IMPACTS')
dataPath = pathlib.Path(BASE_DIR + dataDir + '/ONLY_IMPACTS')
pickles = list(dataPath.glob('*_grs1.pickle'))

modelDir = BASE_DIR + dataDir + '/models'
usePtot = True  # Old boolean for ignoring momentum
pop_names = ['JFC', 'HTC', 'AST', 'OCC', 'Uniform']
populations = []
# Takes about 30 seconds to read all populations in
for p in pop_names:
	print(p)
	populations.append(pop(modelDir = modelDir, pop_type = p, usePtot = True))

print("Reading through pickle files")

logFile = open("readCatalog.log", "w+")


i = 0
# run through the pickles
for p in pickles:
	print(p)

	# identify segment
	segment = str(p.stem[0:10])

	# make plot directory
	plotDir = BASE_DIR + '/plots/' + str(segment)
	if os.path.exists(plotDir):
		pass
	else:
		os.makedirs(plotDir)
	i += 1
	if not os.path.exists(plotDir):
		os.makedirs(plotDir)

	# load GRS1 data and initalize class instance
	chainFile = str(dataPath) + '/' + str(segment) +'_grs1.pickle'
	impact1 = impactClass(chainFile, GRS_num = 1)

	logFile.write(segment + '\n')

	# load GRS2 data
	chainFile = str(dataPath) + '/' + str(segment) +'_grs2.pickle'
	if os.path.exists(chainFile) :
		fid = open(chainFile,'rb')
		impact2 = impactClass(chainFile, GRS_num = 2)
		logFile.write("\t Loaded GRS2" + '\n')

		impact_list = [impact1, impact2]
	else:
		impact_list = [impact1]

	try:
		# Make Dual Corner
		hf = impact1.dualCorner(impact2)
		hf.savefig(plotDir + '/dualCorner.png', format = 'png')
		plt.close(hf)
		print('finished making dual corner plot')

		# Goes through GRS1 and GRS2 
		for impact in impact_list:

			# Skymaps
			# make skymaps for GRS1 and 2
			frames = ['sun', 'sc']
			for frame in frames:
				fig = impact.makeMollweide(frame = frame)
				fig.savefig(plotDir + '/sky_%s%i.png'%(frame, impact.grs), format = 'png')
				plt.close(fig)

			# Flat LPF
			fig = impact.makeFlatLPF(N = 50, scale = 'log')
			fig.savefig(plotDir + '/flat_LPF_log_GRS%i.png' % (impact.grs),format = 'png')
			plt.close(fig)

			impact = impact.SCtoSun()
			fig = impact.plot_populations(populations, norm = True, scale = 'lin', show_impact = True)
			fig.savefig(plotDir + '/pops_grs%i.png'%(impact.grs), format = 'png')
			plt.close(fig)

			if args.gif:
				# 3D LPF
				# Fill directory with gifs
				impact.fillGifDir(BASE_DIR + '/plots/' + str(segment))
				# Make gifs from directory
				impact.gif_maker(BASE_DIR + '/plots/' + str(segment))

				print("Finished Making Gifs")
			logFile.write("All done with segment \n\n")
			
	except Exception as e:
		print("Something Went wrong!!")
		print(e)
		traceback.print_exc()
		logFile.write("Something went wrong in segment " + segment + "\n")

		continue
	
	print("all files succesfully created \n")
print("Done with all!")
logFile.write("Done with all")
logFile.close()
