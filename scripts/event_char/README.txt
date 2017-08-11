nicolepagane|README.txt

#######################
STRUCTURE OF README
1. neccessary packages/software
2. general execution of files **MOST IMPORTANT TO READ IF ANYTHING** 
3. expected mcmc output structure 
4. python scripts
5. file structure

#######################
NECESSARY PACKAGES/SOFTWARE

python 3.6 (but 2.7 also works i think)
anaconda 4.4.0 (should cover most of the necessary modules like the following)
	numpy 1.12.1
	matplotlib 2.0.2
		Basemap 1.1.0
	pandas 0.20.1
corner 2.0.2

#######################
GENERAL EXECUTION OF FILES

each python script is a module that contains functions and data when imported 
if executed at the command line via 'python *script*.py', the script asks the user whether certain functions should be implemented or not
	suggestions are printed to best advise the user how to run the files
	since there are text-prompts during command-line execution (and thus intended to be self-explanatory), it is not entirely necessary to read this whole README
		sparknotes version:
			read the rest of this section, the expected MCMC structure, and the DEPENDENCIES for each python script

if in a python environment, type 'import *script*' and the functions and data are now in your working environment to use as you please 
	type 'help(*script*)' in python for a better idea about which functions and data are available in this module
	if you do prefer to use the specific commands in python rather than executing the file at the command line, i would suggest looking at the bottom of the file (after "if __name__ == '__main__':" to see how/what i intended the functions to be used)

#######################
EXPECTED MCMC OUTPUT STRUCTURE	

The MCMC outputs of micrometeroite impacts on LPF contain the following parameters:
1. Log Likelihood
2. Number of impacts
3. impact time from the start of the segment in seconds (segment start is in the filename)
4. impact momentum in N-s
5. not sure
6. not sure
7. cosine latitude (in spacecraft coordinates)
8. longitude (in spacecraft coordinates)
9. impact face (not sure of the numbering)
10. impact location in x (m)
11. impact location in y (m)
12. impact location in z (m)

this is the assumption of all the scripts so if this order is to be adjusted (as i think tyson said it might to include SNR), i'll have to fix this (namely bc i read in the files differently as the summer progressed as i learned more about python) so it won't be a very simple fix

#######################
PYTHON SCRIPTS

####
catalog_tool.py 
This module contains funtions to read, process, and plot the MCMC impact data. In addition, it develops the html page to visualize the plots for each impact event and appends the catalog. 
it's not really necessary to make plots and webpages when executing due to Sophie's event viewer scripts being much better 
you can run this script to add to the catalog whenever an additional event is found and added to the data directory (since it does not append data that is already existent in the catalog)
DEPENDENCIES:
	data/
		run_e_*GPStime*/ (impact subdirectories should be on tsankawi around local/data/ltpda)
			of which contain several files, but specifically an impactchain.dat file

GENERATES:
	html/ (this is only generated if the user says so when executing the file)
		*GPStime*.html webpages that contain plots of the posterior parameter distributions
		figs/
			*GPStime* (impact subdirectories, from where the *GPStime*.html webpages get the plots)
				jpgs of the posterior distributions 
	catalog.dat (impact catalog with 90% CI parameters in the form [min mid max], where the first 2 lines are headers/comments and the delimiter is a single space) 
	catalog50.dat (impact catalog of parameters at 50% confidence)	
	all_events.html (only if user says so) ~ master webpage from where event webpages w plots can be accessed
	
HOW TO RUN:
	you can just run it at the command line via 'python catalog_tool.py'
	it will prompt the user whether or not they want to generate plots and webpages (you shouldnt since sophies are better)

####
rotation_tool.py
This script contains functions to convert and test for rotations from the LPF SC frame into the meteor frame.
When executed, it prompts the user to either check the z-location (sun) before rotations or to go ahead and perform rotations for all the available MCMC data.

DEPENDENCIES:
	data/
		run_e_*GPStime* (impact subdirectories)
                        of which contain several files, but specifically an impactchain.dat file
	sky_angles/allQuats.txt (mission attitude quaternion files found on tsankawi at /Home/eud/jslutsky/Public/allQuats.mat)

GENERATES:
	sky_angles/meteor/orientation/
		run_e_*GPStime*.dat (meteor frame sky angle data [lat, lon] both in radians w delim = ' ')

HOW TO RUN:
	python rotation_tool.py

####
model_tool.py
This file contains functions to make model inferences on the micrometeorite populations (JFC = Jupiter-Family Comets, HTC = Halley-Type Comets), namely
rates from extrapolated model fluxes given sky location and momentum.
There is a text-prompt function of the rate interpolation and a rate function to use for the MCMC data 
THE MCMC RATE FUNCTION CODE NEEDS MORE WORK BC THE INTERPOLATION IS TOO SLOW 
I NEED TO PROBABLY JUST USE BINARY SEARCH BUT IDK WHY IM BEING STUBBORN AND TRYING OTHER THINGS FIRST

DEPENDENCIES:
        data/
                run_e_*GPStime*/ (impact subdirectories)
                        of which contain several files, but specifically an impactchain.dat file
        sky_angles/
		allQuats.txt (mission attitude quaternion files found on tsankawi at /Home/eud/jslutsky/Public/allQuats.mat)
		meteor/orientation/
			run_e_*GPStime*.dat (meteor frame sky angle data [lat, lon] both in radians w delim = ' ')
	models/
		JFC_30um (JFC model file from Petr and Diego)
		HTC_30um (HTC model file from Petr and Diego)
	
GENERATES:
	models/rates/
		JFC/
			run_e_*GPStime*.dat (expected rates given sky angle and momentum and squared distance result from the interpolation to see how good of a selection the nearest neighbor joining was) delim = ' '	
		HTC/
			run_e_*GPStime*.dat (^ same)
		*** DONT RUN AS IS YET BC IT IS VERY SLOW AND I CAN IMPROVE IT BUT I RAN OUT OF TIME AT THE MOMENT

HOW TO RUN:
	python model_tool.py
		answer prompts (and for now just use the prompt-rate function instead of the MCMC function so that you dont sit around and let the code run forever, sorry again)

####
plot_tool.py
This file contains functions to generate skymaps of 2D histograms of the impact sky locations and some other random plotting tools.
it asks the user whether they want to make individual plots for each event, successive histograms for all the data in order of increasing uncertainty in localization, and what weighting to use
*Once i get the interpolation to run quicker, i'll add another function to make a corner plot to look at correlations between rate, momemtum, and sky angles

DEPENDENCIES:
        data/
                run_e_*GPStime*/ (impact subdirectories)
                        of which contain several files, but specifically an impactchain.dat file
        sky_angles/
                meteor/orientation/
                        run_e_*GPStime*.dat (meteor frame sky angle data [lat, lon] both in radians w delim = ' ')
        *models/rates/
		JFC/
			run_e_*GPStime*.dat
		HTC/
			run_e_*GPStime*.dat        **not yet but when i write the additional corner plot function after more efficient interpolation

GENERATES:
	sky_angles/meteor/figs/
		2d histogram of sky angles (meteor frame) jpgs for each individual event
		series/
			appending successive histograms increasing in location error 
	sky_angles/SC/figs/
		2d histogram of sky angles (SC frame) jpgs for each individual event
                series/
                        appending successive histograms increasing in location error 
	meteorskymap.jpg (all MCMC sky angles in 2d histogram ~ meteor frame)
	SCskymap.jpg (all MCMC sky angles in 2d hist ~ SC)

HOW TO RUN:
	python plot_tool.py
		just answer the prompts

#######################
FILE STRUCTURE 

event_char/

	data/
		run_e_*GPStime*/

	html/
		figs/
			*GPStime*/

	models/	
		rates/
			JFC/
			HTC/

	sky_angles/
		meteor/	
			figs/
				series/
			orientation/
		SC/	
			figs/
				series/
