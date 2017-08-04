Scripts written by Sophie Hourihane (University of Michigan) 
For fufillment of the 2017 UAH/NASA Marshall Heliophysics REU

Currently Located at, Need a LIGO ID 
https://ldas-jobs.ligo-wa.caltech.edu/~sophie.hourihane/REU/

--------------------------------------------o-----
-------------------------------------------o------
------------------------------------------o-------
-------------- ___________________-------o--------
----------- /                      \----o---------
----------/                          \-o----------
------- /                              \-o--------
-------|        LISA PATHFINDER         |--o------
-------|     MICROMETEOROID IMPACTS     |----o----
-------|                                |------o--
-------|                                |--------o
-------|                                |--------- o
--------\                              /----------   o
----------\                          /------------
------------\______________________/--------------                 
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
  

________ INDEX _________

I) Necessary to Read
   	A) Necessary packages or updates
   	B) Changes to scripts before running
	C) Probable errors

II) Script explanations
   Naming conventions and shorthand

   A) makegraphs.py
        - input
		1) Executables  
	- output       
		1)  Flattened LISA Pathfinder
		2)  Momentum Histogram
		3)  Sky Location of Impact
		4)  Likelihood Chain
		5)  Impact Timeline
		6)  CDFs
		7)  Impact Reconstruction
		8)  Noise Constrain
		9)  Noise Constraint Confidence
		10) Gifs

   B) makehtml.py 
   C) lpf_graph_functions.py
   D) patch3dfunc.py
   E) colormap.py
   F) gif_func.py

III) Directory organization



###########################################################################################

--------------------------          I           ----------------------------------------
--------------------------   Necessary to Read  ----------------------------------------

###########################################################################################

__________________________________________________________________________________
                                     - A -
 -------------------------- List of Necessary Packages --------------------------
__________________________________________________________________________________

To run these scripts, you will absolutely need 

numpy (1.13.1)        -The scripts will CRASH with older versions of numpy

the following language and packages may not be totally necessary to update, but these worked

python(2.7.5)

- pandas(0.20.3)
- matplotlib(2.0.2)
- scipy(0.19.1)
- imageio(2.2.0)      -Used for making gifs

If you can't update these because you are on a remote computer, I made a virtual environment.
This documentation was helpful. 
http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/

__________________________________________________________________________________
                                     - B -
 ------------------------ Changes to script before running ----------------------
__________________________________________________________________________________

_______________ makegraphs.py ______________________

inside /REU/scripts/makegraphs.py line 40, marked clearly with ###

1) define where the REU folder is

ex) homedir = '/home/sophie.hourihane/public_html/REU'

-hint) do not put a '/' at the end of 'REU'


_______________ makehtml.py ______________________

inside /REU/makehtml.py line 10, marked clearly with ###

1) define where the REU folder is

ex) homedir = '/home/sophie.hourihane/public_html/REU'

-hint) do not put a '/' at the end of 'REU'


2) define where your website is! 

ex) webaddress = "https://ldas-jobs.ligo-wa.caltech.edu/~sophie.hourihane/REU"  

-hint) do not put a '/' at the end of 'REU'


_____________________________________________________________________
                         - C -               
------------- Errors That Wouldn't Suprise Me ---------------------
_____________________________________________________________________

... There may be some whitespace issues and I apologize in advance, I was working from 
2 different computers and they dissagreed about what a tab was. This might be yikes to fix. 

These progams will throw back some errors, but most of these have to do with me passing
empty arrays, arrays with a single value, or masked arrays. This is also why you need the newest version of numpy :) 
If the program doesn't crash, the error is just a warning that it hates what its doing.

Some of the fontsizes may look bad. If you want it changed globally throw this at the beginning of makegraphs.py after import statements
	ex)
	# Set Global fontsize 
        matplotlib.rcParams.update({'font.size': 12})
	
Except you can put anything instead of 12
 


##########################################################################################

---------------------------------        II           ---------------------------------
--------------------------------- Script explanations ---------------------------------

##########################################################################################

______________________________ Naming Conventions __________________________

In script names
- if the script ends with func.py, then it is called by makegraphs.py, and need not be changed
- makegraphs.py and makehtml.py are the only functions that are called individually

Inside of scripts
... messy to say the least, but here are some that are pretty consistent

 currun         = impact GPS time, 10 digits, (curr meaning current, bad convention)
 mednmom        = median momentum 
 LPF or lpf     = LISA pathfinder
 df_(anything)  = a pandas dataframe, always filled with LPF data

_Apologies for bad naming_

confidence really means credible

_________________________________      A      ______________________________
_________________________________makegraphs.py______________________________

How to run it:

	call python
		ex) python REU/scripts/makegraphs.py

What it Does:
	Without any extra commands this will make the most general plots,will run for ALL times, and will probs take 1/2 hour to hour to run
		Plots)  - Impact Timeline
			- Momentum and Fractional Error CDF
			- Momentum Histograms  	
			- Sky Location of Impact
			- MCMC Log Likelihood Chain
			- Model Reconstruction in Acceleration, amplitude vs time
			- Impact map of the LPF, Flattened 
Different Executables:

	To run one of these executables, type
		python REU/scripts/makegraphs.py (EXECUTABLE COMMAND(s))
	--------------------------------------------------------------
	
	a) Help Menu
		
		This will show you a help menu! 

		COMMAND) -h OR --help
		
		Examples)
			python REU/scripts/makegraphs.py --help
			python REU/scripts/makegraphs.py -h
			
	
	b) Run only on a single GPS time

		COMMAND) -d GPS_TIME, --directory GPS_TIME
		
		Examples) 
			Produce default graphs for single directory 
				python REU/scripts/makegraphs.py --directory 1144228507
				python REU/scripts/makegraphs.py -d 1144228507
					
			You can also add this command on top of others:

			Will produce a gif for the GPS time 1144228507
				python REU/scripts/makegraphs.py -d 1144228507 --gif
	
	c) Produce all possible graphs (including gifs)
		
		COMMAND) -e, --everything
		
		Examples)
			Make every. single. graph. (Will take all work day) 
				python REU/scripts/makegraphs.py -e

			Make all the graphs for a single GPS Time (Warning, Cumulative Plots will be overwritten incorrectly)	
				python REU/scripts/makegraphs.py -d 1144228507 -e
	
	d) Make the Noise Parameterizing Plots
		
		COMMAND) --noise

		Examples)
			Makes the plots parameterizing the noise of all GPS times
				python REU/scripts/makegraphs.py --noise

			Make all the noise graphs for a single GPS Time (Warning, Cumulative Plots will be overwritten incorrectly)	
				python REU/scripts/makegraphs.py -d 1144228507 --noise


	e) Make the 3D LPF Gifs
		
		COMMAND) -g, --gifs, --gif

		Examples)
			Makes every gif, Will run all work day. 
				python REU/scripts/makegraphs.py --gif
			
			Make a gif for a specified GPS time
				python REU/scripts/makegraphs.py -d 1144228507 --gif

	f) Run the cumulative plots
		
		COMMAND) -cumul, --cumulative

		Examples)
			Makes only the timeline, momentum histograms, and the momentum confidence intervals
				python REU/scripts/makegraphs.py --gif
			

	

                 __________________ Outputs ____________________

1) Flattened LISA Pathfinder 

	- filename
		lpf_flat_( **  )_currun_mednmom.png 

	- SWITCH
		flatten = True

	- What is it? 
		The flattened version of LISA pathfinder to cut up on paper
		A 2d Histogram of all faces of the LPF, rotated to show continuity between faces
		Legend refers to the percent of impacts on that "face" of the lpf

	- Function Location
		located in REU/scripts/lpf_graph_functions.py
		function is titled : flatten_LPF

	- File Location
		REU/runs/run_e_currun/images_currun_mednmom/

	- ** 
		lin medning linear scale in histogram, 
		log medning log scale in histogram

	- Possible Future Problems
		 None at the moment

	- Possible Future Improvements
		The length of each bin is close, but not exactly the same, fixing this would be nice


2) Momentum Histogram

	-  filename
		momentum_currun_mednmom.png

	- SWITCH 
		momentum_hist = True

	- What is it?
		histogram of MCMC genated momenta
		number of bins decided with sqrt(length(momenta))

	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/run_e_currun/images_currun_mednmom/

	- Possible Future Problems
		Some weird opacity issues. I think its a matplotlib problem

	- Possible Future Improvements
		A smarter way to make bins, right now sqrt(length)  


3) Sky Location of Impact

	-filename
		skyloc_currun_mednmom.png

	- SWITCH
		skyloc = True

	- What is it?
		makes a 2d Histogram of the origin (sky location) of the impact in spacecraft coordinates
		x = longitude
		y = cosine of latitude 

	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/run_e_currun/images_currun_mednmom/

	- Possible Future Problems
		None at the moment

	- Possible Future Improvements
		A smarter way to make bins, currently sqrt(length(longitude)) 


4) MCMC Log Likelihood Chain

	-filename
		logprob_chain_currun_mednmom.png

	- SWITCH
		show_path = True

	- What is it?
		makes a 2d Histogram of the sky location in spacecraft coordinates
		x = longitude
		y = cosine of latitude 

	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/run_e_currun/images_currun_mednmom/

	- Possible Future Problems
		None at the moment

	- Possible Future Improvements
		None at the moment


5) Impact Timeline

	-  filename
		timeline.png

	- SWITCH 
		timeline = True

	- What is it?
		A timeline of the impacts, UTC Time as well as sections that we have looked at. 

	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/total_hists/

	- Possible Future Problems
		The catalogue filenames (GPS and length of time sections that we have looked at)
		 are hardcoded based on splicing strings. 
		If the looking time changes past 4 digits, this will no longer be valid. Sorry!  

	- Possible Future Improvements
		A way around the above problem. 

6) CDFs (Momentum Credible and Fractional error) 

	-  filenames
		mom_90_confidence.png                Fractional Error CDF
		mom_confidence.png                   90% 50% Median CDF for momentum 

	- SWITCH 
		mom_confidence = True

	- What is it?
		CDFs for the momentum and for the Fractional error 
	
	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/total_hists/

	- Possible Future Problems
		None at the moment

	- Possible Future Improvements
		None at the moment 


7) Impact Reconstruction 

	-  filenames
		impact_pulse_GPS_Momentum.png

	- SWITCH 
		show_impacts = True

	- What is it?
		Shows Model and data agreement in the reconstruction of the impact
		 for each Degree of freedom 
 
	- Function Location
		located in makegraphs.py

	- File Location
		REU/runs/run_e_GPS/images_GPS_Momentum/

	- Possible Future Problems
		None at the moment

	- Possible Future Improvements
		None at the moment 

8) Noise Constraint  

	-  filenames
		amp_hist.png
		params_hist.png
		q_hist.png	

	- SWITCH 
		make_peak   = True 
		make_amp    = True
		make_q      = True
		
	- What is it?
		make_peak - Makes a plot of where the peaks in the noise are for each DOF 
		make_amp  - Makes a 2d Histogram of where the peaks in the Freqs vs amps are
		make_q   - Makes a 2d Histogram of where the peaks in the Freqs vs q are
	
	- Function Location
		located in freqlines.py
		called by makegraphs.py

	- File Location
		REU/runs/run_e_GPS/images_GPS_Momentum/hists/

	- Possible Future Problems
		These functions are not refined well, and its unclear how well they worked at the beginning.  

	- Possible Future Improvements
		Make these work at all? They are still included but its hard because the time analyzed 
		was small which makes the lines bad

9) Noise Constraint Confidence

	-  filenames
		dofsnoise_currunGPS_mednmomMOMENTUM.png     .... I know its a bad name
	
	- SWITCH 
		make_flines = True
	
	- What is it?
		Makes 90% credible intervals around the noise as well as the noise 
	
	- Function Location
		Located in lpf_graph_functions.py under mean_confidence_interval
		called by makegraphs.py

	- File Location
		REU/runs/run_e_GPS/images_GPS_Momentum/

	- Possible Future Problems
		Hard coded to skip first 50 of 100 files. If more / less are added, could be a problem because this will change

	- Possible Future Improvements
		Fix above

10) Gifs

	-  filenames
		pf_log3d_ANGLE.png                      -- There are many of these, each one is a frame in the gif 
		gif_runGPS_momMOMENTUM_(top or bot).gif -- This is the actual gif file  
	
	- SWITCH 
		makegifs
		make_gifs_all
	
	- What is it?
		makegifs      makes and fills the image directories with the images for the gifs, this takes all day
		make_gifs_all makes the actual .gif files	
	
	- Function Location
		makegifs - patch3dfunc.py
		make_gifs_all - gif_fun.py

	- File Location
		makegifs - REU/runs/run_e_GPS/images_GPS_Momentum/gifs/log_3dlpf_momMomentum_runGPS_(bot or top)/
		make_gifs_all -REU/runs/run_e_GPS/images_GPS_Momentum/
	
	- Possible Future Problems
		None at the moment
	- Possible Future Improvements
		Make this faster!! 

ERRORS TO IGNORE :
        (More than just this LOL)
 	
	lib/python2.7/site-packages/matplotlib/colors.py:469: UserWarning: Warning: converting a masked element to nan.
	  xa = np.array([X])

_________________________________      B        _____________________________
_________________________________  makehtml.py  _____________________________

How to run it:
	python makehtml.py

What it Does:
	Makes the html directories and all the html pages
	
Before you run: 
	make sure that you change the webaddress and the homedir ! 	

____________________________           C            _________________________
____________________________ lpf_graph_functions.py _________________________


How to run it:
	This function is called automatically by makegraphs.py

What it Does:
	This has the more complicated functions that took up too much space 
	to be included in makegraphs.py

Functions included:
	Flatten LPF (and associated functions)
	Some old ones that are no longer in use (createLPF)
	90% Credible interval for Noise (... titled ploterror)
	


____________________________           D            _________________________
____________________________     patch3dfunc.py     _________________________


How to run it:
	This function is called automatically by makegraphs.py

What it Does:
	Has the function that creates the 3d LPF, a long function

Functions included:
	3D LPF and associated functions


____________________________           E            _________________________
____________________________      colormap.py       _________________________


How to run it:
	This function is called automatically by makegraphs.py

What it Does:
	Just defines the parula colormap. It's not a part of matplotlib :( 	


_________________________________      F      ______________________________
_________________________________ gif_func.py ______________________________


How to run it:
	This function is called automatically by makegraphs.py

What it Does:
	Uses imageio (a library)	
	Makes 2 gifs per impact, 15 degrees above and below the x axis
	



######################################################################################

---------------------------------- Directory Organization ---------------------------

######################################################################################

This is how these files and directories should be organized (Most of these are made within the scripts)

  ______________LEGEND ____________

  *name*      generated by running the scripts
  :           list of contents and or subdir
  /           directory
  _________________________________

                                                    /REU
          __________________________________________ : _______________________________________
           :                          :                            :                            
        /scripts                     /runs                    */html scripts*               
           :		              :                            :                            
 graph making scripts                 :               *.html files for webpage*       
	                              :
	______________________________:____________
        :		   :                      :
	:		   :                      :
catalogfiles.txt   /run_e_(gpstime)         */total_graphs*		      
			 :                          :
		_________:________                  :__________________________
               :                  :                                           :
  .dat files from LPF             :                                 *graphs using data from all runs*
                         */images_(gps time)_(medn momentum)*
                                  :
            _____________________ : _________________
           :                                         :
 *.png and .gifs generated for this impact*       */gifs*
                                           __________:________
					 :                    :
				   */..top*                  */..bot*
                                      :                         :
		   .png view from 15 deg above 3D LPF    .png viewed from 15 deg below 3D LPF  

					








