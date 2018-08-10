import pandas as pd
import numpy as np
import pathlib
import os

"""
This function takes Petr's population models from
being cumulative to being bins
Sophie Hourihane 7/2018
"""
p = pathlib.PurePath(os.getcwd())
BASE_DIR = str(p.parent)
dataPath = pathlib.Path(BASE_DIR + '/data')

dataDir = dataPath
pops = ['AST', 'JFC', 'HTC', 'OCC']
for pop in pops:
	print(pop)
	if pop == 'AST':
		dataFile = 'population_models/AST_NEW_impulse.dat'
	else:
		dataFile = 'population_models/' + pop + '_impulse.dat'


	Ptots = [1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
	df = pd.read_csv(str(dataDir) + '/' + dataFile, 
					 names = ['lon', 'lat', 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3], sep = '\s+')

	df_new = df.copy()


	print('Before \n', df_new.loc[[8]])
	for j in range(len(Ptots) - 1):
		df_new[Ptots[j]] = df[Ptots[j]] - df[Ptots[j + 1]]
	print('After \n', df_new.loc[[8]])

	print('should be same number')
	print(sum(df[1e-7].values))
	print(sum(df_new[Ptots].sum(axis = 1).values))

	print('\n')

	new_file = open(str(dataPath) + '/models/' + pop + '_impulse.dat', 'w+')
	for i in range(len(df_new['lon'])):
		for j, p in enumerate(Ptots):
			new_file.write('%.f\t%.f\t%1.1e\t%1.6e\n'%(df_new['lon'][i], df_new['lat'][i], p, df_new[p][i]))
