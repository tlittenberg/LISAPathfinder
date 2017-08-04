import numpy as np
import glob
import os
from os import listdir
from os.path import isfile, join
from datetime import datetime, timedelta
import re

########################################### Change These Lines ! ###############################################
homedir      = '/home/sophie.hourihane/public_html/REU'                       #Absolute Location of REU folder
webaddress = "https://ldas-jobs.ligo-wa.caltech.edu/~sophie.hourihane/REU"    #Web Address

#Do not put / at the end of these lines or stuff won't run 
################################################################################################################

# Sets the style for the whole website
style = """<!DOCTYPE html>
        <html>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <link rel="stylesheet" href="https://www.w3schools.com/w3css/4/w3.css">
        <head>
        <style>

        body {
          font-family: "Lato", sans-serif;
          margin:0;
          background-color:#FCF4D9 }

        .topnav {
          font-family: "Lato", sans-serif;
          overflow: hidden;
          background-color: #555;
         
        }

        div.gallery.noise {
            margin: 13px;
            border: 1px solid #ccc;
            float: left;
            width: 600px;
            height: auto;
        }

        div.gallery.small {
            margin: 5px;
            border: 1px solid #ccc;
            float: left;
            width: 600px;
            height: auto;
        }
        div.gallery:hover {
            border: 1px solid #777;
        }

        div.gallery img {
            width: 100%;
            height: auto;
        }


        .topnav a {
          font-family: "Lato", sans-serif;
          float: left;
          display: block;
          color: #f2f2f2;
          text-align: center;
          padding: 14px 16px;
          text-decoration: none;
          font-size: 17px;
        }

        .topnav a:hover {
          font-family: "Lato", sans-serif;
          background-color: #ddd;
          color: black;
        }

        .topnav a.active {
                font-family: "Lato", sans-serif;
            background-color: #D55B3E;
            color: white;
        }

    * {
      box-sizing: border-box;
      font-family: "Lato", sans-serif;
    }

    #myInput {
      background-position: 10px 12px;
      background-repeat: no-repeat;
      width: 75%;
      font-family: "Lato", sans-serif;
      font-size: 16px;
      padding: 12px 20px 12px 40px;
      border: 1px solid #ddd;
      margin-bottom: 12px;
    }

    #myUL {
      list-style-type: none; 
      font-family: "Lato", sans-serif;
      width: 75%;
      padding: 0;
      margin: 0;
    }

    #myUL li a {
      border: 1px solid #ddd;
      margin-top: -1px; /* Prevent double borders */
      background-color: #f6f6f6;
      padding: 12px;
      text-decoration: none; 
      font-size: 18px;
      font-family: "Lato", sans-serif;
      color: black;
      display: block
    }

    #myUL li a.header {
      background-color: #e2e2e2;
      cursor: default;
    }

    #myUL li a:hover:not(.header) {
      background-color: #eee;
    }

        a {
            text-decoration: none;
            display: inline-block;
            padding: 8px 16px;
        }

        a:hover {
            background-color: #ddd;
            color: black;
        }

        .previous {
            background-color: #f2f2f2;
            color: black;
        }

        .next {
            background-color: #D55B3E;
            color: white;
        }

        .round {
            border-radius: 50%;
        }
        
	.dropdown {
	    float: left;
	    overflow: hidden;
	}

	.dropdown .dropbtn {
	    font-size: 16px;    
	    border: none;
	    outline: none;
	    color: white;
	    padding: 14px 16px;
	    background-color: inherit;
	}

	.container a:hover, .dropdown:hover .dropbtn {
	    background-color:#D55B3E;
	}

	.dropdown-content {
	    display: none;
	    position: absolute;
	    background-color: #f9f9f9;
	    min-width: 160px;
	    box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
	    z-index: 1;
	}

	.dropdown-content a {
	    float: none;
	    color: black;
	    padding: 12px 16px;
	    text-decoration: none;
	    display: block;
	    text-align: left;
	}

	.dropdown-content a:hover {
	    background-color: #ddd;
	}

	.dropdown:hover .dropdown-content {
	    display: block;
	}

        </style>
        </head>"""

topnav =  """<body>

        <div class="topnav">
          <a class="active" href="index.html">Home</a> 
          <a href="momenta.html">Momenta Search</a>
	  <a href="time.html">Time Search</a>
	<div class="dropdown">
		<button class="dropbtn">Time Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_runs.html">Gifs</a>
			<a href="momenta_graphs_runs.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>	
			<a href="flat_LPF_runs.html">Flat LPFs</a>
			<a href="likelihood_runs.html">Likelihood</a>
			<a href="skyloc_runs.html">Sky Location</a>
			<a href="impact_runs.html">Impact Model</a>
		</div></div>
	<div class="dropdown">
		<button class="dropbtn">Momenta Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_moms.html">Gifs</a>
			<a href="momenta_graphs_moms.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>
			<a href="flat_LPF_moms.html">Flat LPFs</a>
			<a href="likelihood_moms.html">Likelihood</a>
			<a href="skyloc_moms.html">Sky Location</a>
			<a href="impact_moms.html">Impact Model</a>
		</div></div>
	</div> """


#sorts List
def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

# Converts GPS time to UTC
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




# Sorts List
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    s = text 
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    final_list = [float(x) for x in re.findall(match_number, s)]
    
    return final_list


def get_the_subdir(a_dir):
    subdir = []
    names  = []
    #print('directory = %s'%(a_dir))
    for name in os.listdir(a_dir):
        #print('name of directory = %s'%(name))
        if os.path.isdir((os.path.join(a_dir,name))):
            names.append(name)
            subdir.append((os.path.join(a_dir, name)))
    subdir.sort(key=natural_keys)
    names.sort(key=natural_keys)
    return subdir, names


def makeindex(webaddress,title, scriptsdir):
    filename = ('%s/%s.html'%(scriptsdir,title))
    #filename = ('%swebsite/public_html/html_scripts/%s.html')%(homedir,title)
    f = open(filename,'w+') 

    header = style 

    test = """<body>

        <div class="topnav">
          <a class="active" href="index.html">Home</a> 
          <a href="momenta.html">Momenta Search</a>
	  <a href="time.html">Time Search</a>
	<div class="dropdown">
		<button class="dropbtn">Time Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_runs.html">Gifs</a>
			<a href="momenta_graphs_runs.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>	
			<a href="flat_LPF_runs.html">Flat LPFs</a>
			<a href="likelihood_runs.html">Likelihood</a>
			<a href="skyloc_runs.html">Sky Location</a>
			<a href="impact_runs.html">Impact Model</a>
		</div></div>
	<div class="dropdown">
		<button class="dropbtn">Momenta Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_moms.html">Gifs</a>
			<a href="momenta_graphs_moms.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>
			<a href="flat_LPF_moms.html">Flat LPFs</a>
			<a href="likelihood_moms.html">Likelihood</a>
			<a href="skyloc_moms.html">Sky Location</a>
			<a href="impact_moms.html">Impact Model</a>
		</div></div>
	</div>"""
    navigation = topnav +"""</div>
        <div style="padding-left:16px">
    		<h2> LISA Pathfinder Micrometeoroid Impacts </h2>
	<p> A spacecraft designed in Europe <br>
	  Whose data impacts interrupt <br>
	  Till Thorpe saw one day <br>
	  It need not be that way. <br>
	  We just need that Littenberg up </p>
    
	<p> This mission created by ESA <br>
	  Proved tech was ready for LISA  <br>
	  Impulse sensitivity <br>
	  Gave us a proclivity <br>
	  To research what comets released us </p>


    </body>
    </html>"""
    whole = header + navigation 

    f.write(whole)
    f.close()

def makeimgpage(webaddress,title,filenames,i, img_directory, scriptsdir):
	##print(title[i])
        #print('making directory %s' %(directory))
		
        filename = ('%s/%s.html'%(scriptsdir,filenames[i]))
	#filename = ('%s/website/public_html/html_scripts/%s.html')%(homedir,title[i])
        f = open(filename,'w+')

	wraptitle = style 
        
        nexti = i+1
        previ = i-1
        if previ == 0:
            previ = i
        if nexti == len(filenames):
            nexti = i
            #nexti = len(title)-2
        filenamenext = ('%s/html_scripts/%s.html'%(webaddress,filenames[nexti]))
        filenameprev = ('%s/html_scripts/%s.html'%(webaddress,filenames[previ]))
	#filenamenext = ('%s/website/public_html/html_scripts/%s.html')%(webaddress,title[nexti]) 
        #filenameprev = ('%s/website/public_html/html_scripts/%s.html')%(webaddress,title[previ])
        body = """<body>

        <div class="topnav">
          <a class="active" href="index.html">Home</a> 
          <a href="momenta.html">Momenta Search</a>
	  <a href="time.html">Time Search</a>
	<div class="dropdown">
		<button class="dropbtn">Time Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_runs.html">Gifs</a>
			<a href="momenta_graphs_runs.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>	
			<a href="flat_LPF_runs.html">Flat LPFs</a>
			<a href="likelihood_runs.html">Likelihood</a>
			<a href="skyloc_runs.html">Sky Location</a>
			<a href="impact_runs.html">Impact Model</a>
		</div></div>
	<div class="dropdown">
		<button class="dropbtn">Momenta Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_moms.html">Gifs</a>
			<a href="momenta_graphs_moms.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>
			<a href="flat_LPF_moms.html">Flat LPFs</a>
			<a href="likelihood_moms.html">Likelihood</a>
			<a href="skyloc_moms.html">Sky Location</a>
			<a href="impact_moms.html">Impact Model</a>
		</div></div>
	</div>
        <div style="padding-left:16px">
          <h2>%s</h2>
        

        <a href="%s" class="previous round">&#8249;</a>
        <a href="%s" class="next round">&#8250;</a>
  
        </div>
        """%(title[i], filenameprev,filenamenext)



#       <div class="w3-bar w3-teal">
# """
        #for path, subdirs, files in os.walk(directory[0]):
    #           for name in files:
#                       print(name)
#                       print os.path.join(path,name)   
        
        #"""<a href="#" class="w3-bar-item w3-button w3-mobile">London</a>"""
        images = '' 
        data_noise = []
        mom_img = []
        lpf_flat = []
        logprob = []
        skyloc = []
	gif = []
	hist_list = []
	histfiles = []
	impact = []
        i = 0
        for direc in img_directory: 
	    
            onlyfiles = [fil for fil in listdir('%s'%(direc)) if isfile(join('%s'%(direc),fil))]
	    histdir = direc+'/hists'
	    
	    if os.path.exists(histdir):
          	histfiles = [fil for fil in listdir('%s'%(histdir)) if isfile(join('%s'%(histdir),fil))]
	    #for direc in directory:
	    direc = direc.replace(homedir,'')
	    direc =webaddress + direc 
            
	    histdir = histdir.replace(homedir,'')
	    histdir = webaddress +  histdir
	    for files in histfiles:
	    	#print('#######################################' + direc+'/'+files)
		hist_list.append("""<img src='%s/%s' alt='%s'>"""%(histdir,files, files))    
	    #print('going through files in %s'%(direc))
            #direc = direc.replace(homedir,"")
	    #direc = webaddress+'/'+direc
	    #print(direc) 
            #print(onlyfiles)
            for files in onlyfiles:
                #print(files)
                if 'dofs' in files[:4]:
                    data_noise.append("""<img src='%s/%s' alt='%s'>"""%(direc,files, files))    
#                
                if 'momentum' in files[:8]: 
                    mom_img.append( """<img src='%s/%s' alt='%s'>"""%(direc,files, files))      
#                
                if 'lpf_flat' in files[:8]:  
                    lpf_flat.append( """<img src='%s/%s' alt='%s'>"""%(direc,files, files))     
                
                if 'logprob' in files[:7]: 
                    logprob.append( """<img src='%s/%s' alt='%s'>"""%(direc,files, files))      
		if 'impact' in files[:6]:
                    impact.append( """<img src='%s/%s' alt='%s'>"""%(direc,files, files))      
                if 'skyloc' in files[:6]:
                    skyloc.append("""<p><img src='%s/%s' alt='%s'></p>"""%(direc,files, files)) 
                if 'gif' in files[:3]:
                    gif.append("""<p><img src='%s/%s' alt='%s'></p>"""%(direc,files, files)) 
#                else:
#                    images += """<div class = "gallery"> 
#                                    <img src='%s/%s' alt='%s' width = "600" height ="400">
#                                </div> """%(direc,files, files)        
                    
            i += 1
        
        
        
        smallray = []
        smallray.extend(mom_img)
        smallray.extend(logprob)
        smallray.extend(skyloc)
        small = ''
        
	flat = ''
	for d in sorted(lpf_flat):
            flat += """<div class = "gallery noise"> %s </div>"""%(d)
        images += flat
	#logprob = logprob.sort()
	lprob = ''
	for d in logprob:
            lprob += """<div class = "gallery noise"> %s </div>"""%(d)
	
	images += lprob
	
	sky = ''
        for d in skyloc:
            sky += """<div class = "gallery noise"> %s </div>"""%(d)

	images += sky     
	mom_list = ''
	
	for d in mom_img:	
            mom_list += """<div class = "gallery noise"> %s </div>"""%(d)
	images += mom_list
	impact_list = ''
	for d in impact:	
            impact_list += """<div class = "gallery noise"> %s </div>"""%(d)
	images += impact_list

	giflist =''
	for g in gif:	
            giflist += """<div class = "gallery noise"> %s </div>"""%(g)
        images += giflist
        
	noise_dofs = ''
        for d in data_noise:
            noise_dofs += """<div class = "gallery noise"> %s </div>"""%(d)
            
        images += noise_dofs
        
        hist = '' 
        
	for h in hist_list:
            hist += """<div class = "gallery noise"> %s </div>"""%(h)
	images += hist

	#print(images) 
        #<a href="#" class="w3-bar-item w3-button w3-mobile">Paris</a>
        #<a href="#" class="w3-bar-item w3-button w3-mobile">Tokyo</a>
        
        end = """
        </body>
        </html>"""
   
        whole = wraptitle + body + images + end 
        f.write(whole)
        f.close()
        return filename, giflist, mom_list, flat, hist, sky, lprob, noise_dofs, impact_list

def makegraphpg(homedir,graphtitle,title, graphlist):
#        print('making directory %s' %(directory))
        filename = ('%s/%s.html'%(scriptsdir,graphtitle))
	#filename = ('%s/website/public_html/html_scripts/%s.html')%(homedir,graphtitle)
        f = open(filename,'w+')

        wraptitle = style 
        
	body = """<body>

        <div class="topnav">
          <a class="active" href="index.html">Home</a> 
          <a href="momenta.html">Momenta Search</a>
	  <a href="time.html">Time Search</a>
	<div class="dropdown">
		<button class="dropbtn">Time Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_runs.html">Gifs</a>
			<a href="momenta_graphs_runs.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>	
			<a href="flat_LPF_runs.html">Flat LPFs</a>
			<a href="likelihood_runs.html">Likelihood</a>
			<a href="skyloc_runs.html">Sky Location</a>
			<a href="impact_runs.html">Impact Model</a>
		</div></div>
	<div class="dropdown">
		<button class="dropbtn">Momenta Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_moms.html">Gifs</a>
			<a href="momenta_graphs_moms.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>
			<a href="flat_LPF_moms.html">Flat LPFs</a>
			<a href="likelihood_moms.html">Likelihood</a>
			<a href="skyloc_moms.html">Sky Location</a>
			<a href="impact_moms.html">Impact Model</a>
		</div></div>
	</div>

    	<div style="padding-left:16px">
	<h2> %s </h2>
	</div>
	%s
	
        """%(title,graphlist) 
        end = """
        </body>
        </html>"""
   
        whole = wraptitle + body + end 
        f.write(whole)
        f.close()

def makelistpage(webaddress,filenames, title,search, capital, scriptsdir):
    filename = ('%s/%s.html'%(scriptsdir,title))
    #filename = ('%swebsite/public_html/html_scripts/%s.html')%(scriptsdir,title)
    f = open(filename,'w+')
    
    header = style

    navigation = topnav + """ 

    <div style="padding-left:16px">
    """
    
    searchlist = """<h2>Search by %s</h2>

    <input type="text" id="myInput" onkeyup="myFunction()" placeholder="Search %s.." title="Search %s">
    <ul id="myUL">"""%(capital, capital, capital)


    test = """<body>

        <div class="topnav">
          <a class="active" href="index.html">Home</a> 
          <a href="momenta.html">Momenta Search</a>
	  <a href="time.html">Time Search</a>
	<div class="dropdown">
		<button class="dropbtn">Time Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_runs.html">Gifs</a>
			<a href="momenta_graphs_runs.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>	
			<a href="flat_LPF_runs.html">Flat LPFs</a>
			<a href="likelihood_runs.html">Likelihood</a>
			<a href="skyloc_runs.html">Sky Location</a>
			<a href="impact_runs.html">Impact Model</a>
		</div></div>
	<div class="dropdown">
		<button class="dropbtn">Momenta Graphs</button>	
		<div class="dropdown-content">
			<a href="gifs_moms.html">Gifs</a>
			<a href="momenta_graphs_moms.html">Momenta Hists</a>
			<a href="noise_runs.html">Noise</a>
			<a href="flat_LPF_moms.html">Flat LPFs</a>
			<a href="likelihood_moms.html">Likelihood</a>
			<a href="skyloc_moms.html">Sky Location</a>
			<a href="impact_moms.html">Impact Model</a>
		</div></div>
	</div>

    <div style="padding-left:16px">
    """
    
    searchlist = """<h2>Search by %s</h2>

    <input type="text" id="myInput" onkeyup="myFunction()" placeholder="Search %s.." title="Search %s">
    <ul id="myUL">"""%(capital, capital, capital)
    
    listpart = ""
    index = 0
    for fname in filenames:
	wherefile = fname.replace(homedir+'/',"")	 
        partlink = wherefile.replace("html_scripts/","")
	if 'momenta' in title:
		units = 'kg m/s'
		partlink = partlink.replace('mom_','')
	else:
		units = 's'
		partlink = search[index]#partlink.replace('run_','')
		index += 1
	namelink = partlink.replace('.html','')
	listpart += """<li><a href="%s/%s">%s</a></li>"""%(webaddress,wherefile,namelink+" "+units)

    end = """</ul>
    
    <script>
    function myFunction() {
        var input, filter, ul, li, a, i;
        input = document.getElementById("myInput");
        filter = input.value.toUpperCase();
        ul = document.getElementById("myUL");
        li = ul.getElementsByTagName("li");
        for (i = 0; i < li.length; i++) {
            a = li[i].getElementsByTagName("a")[0];
            if (a.innerHTML.toUpperCase().indexOf(filter) > -1) {
                li[i].style.display = "";
            } else {
                li[i].style.display = "none";

            }
        }
    }
    </script>

    </body>
    </html>"""
    whole = header + navigation + searchlist + listpart + end

    f.write(whole)
    f.close()



directoryruns = homedir + '/runs'                #Location of the run directories
scriptsdir = homedir + '/html_scripts'           #Location of the Scripts directories

if not os.path.exists(scriptsdir):
                os.mkdir(scriptsdir)

makeindex(webaddress,'index', scriptsdir)


#Get all the run_e_ directories
directory_runs, nameruns = get_the_subdir(directoryruns)

imgdir = []		   # Where all the images live, under each run 
mednmoms = []              #List of all the medn momenta
curruns  = []		   #List of all the GPS Times

for d in directory_runs: #There should only be 1 subdirectory, but specifiying that 'images' is in it
	img, name = get_the_subdir(d)
	for i in img:
		if 'images' in i:
			imgdir.append(i)
			mednmoms.append('mom_'+i[-10:])
			curruns.append('run_'+d[-10:])



# Parameterize lists of graphs on momentum sorted pages
i = 0
mom_filenames = []   # Holds all the momentum filenames

gifs = ''
momenta_graphs = ''
flat_lpfs = ''
noise = ''
skyloc = ''
logprob = ''
dofs = ''
impacts = ''
namemoms = mednmoms

#------sort directories by momentum-------#
mednmoms_sort = [] # Sorted by median list
for w in range(len(mednmoms)):
	mednmoms_sort.append(float(mednmoms[w][4:])) # Float the momentum, skip first 4 

mednmoms_sort = np.asarray(mednmoms_sort)
mom_sort = mednmoms_sort.argsort()          # This argument sorts the directories by time

#------sort directories by time-------#
currun_sort = [] #np.zeros_like(curruns)
for w in range(len(curruns)):
	currun_sort.append(curruns[w][4:])
currun_sort = np.asarray(currun_sort)
run_sort = currun_sort.argsort()


#---Makes Graph Pages Sorted by momenta---#
namemoms = np.asarray(namemoms)
imgdir = np.asarray(imgdir)
i = 0

for direc in imgdir[mom_sort]: 
    print(direc,mednmoms[i])
    mom_files, gif_list, mom_list,flat_list,noise_list, sky_list, logprob_list, dofs_list, impact_list = (makeimgpage(webaddress, namemoms[mom_sort], namemoms[mom_sort],i,[direc], scriptsdir))
    mom_filenames.append(mom_files)
    gifs += gif_list
    momenta_graphs += mom_list
    flat_lpfs += flat_list
    noise += noise_list + dofs_list
    skyloc += sky_list
    logprob += logprob_list
    impacts += impact_list
    i += 1

#Makes momentum sorted graph pages
makegraphpg(scriptsdir,'gifs_moms','3D Gifs by Momentum', gifs)
makegraphpg(scriptsdir,'flat_LPF_moms','Impact Location by Momentum', flat_lpfs)
makegraphpg(scriptsdir,'noise_moms', 'Peak Noise Location by Momentum',noise)
makegraphpg(scriptsdir,'momenta_graphs_moms','Momenta Distribution by Momentum', momenta_graphs)
makegraphpg(scriptsdir,'skyloc_moms', 'Sky Location by Momentum',skyloc)
makegraphpg(scriptsdir,'likelihood_moms', 'Likelihood Chain by Momentum',logprob)
makegraphpg(scriptsdir,'impact_moms', 'Impact Model by Momentum',impacts)
#makegraphpg(homedir,'noise_dofs_moms', dofs)


# adds UTC to page
subdirecrun = imgdir
for i in range(len(curruns)):
	realtime = curruns[i]
	realtime = realtime.replace('run_','') #nameruns[i].replace('run_e_','')
	#print(nameruns[i])
	gpstime = int(float(realtime))
	nameruns[i] = "%s, %s"%(gps2utc(0,gpstime),realtime)

i = 0


run_filenames = []
gifs = ''
momenta_graphs = ''
flat_lpfs = ''
noise = ''
skyloc = ''
logprob = ''
dofs = ''
impacts = ''

#print('###########################',curruns, nameruns)
for direc in imgdir:
    run_files, giflist, momlist, flatlist, noiselist, sky_list, logprob_list, dofs_list, impact_list = (makeimgpage(webaddress,nameruns,curruns,i,[direc], scriptsdir))
    run_filenames.append(run_files)
    gifs += giflist
    momenta_graphs += momlist
    flat_lpfs += flatlist
    noise += noiselist + dofs_list
    skyloc += sky_list
    logprob += logprob_list
    impacts += impact_list
#    dofs += dofs_list
    i += 1
#makegraphpg(gifs, 'gifs',
makegraphpg(scriptsdir,'gifs_runs', '3D Impact Location by Time', gifs)
makegraphpg(scriptsdir,'flat_LPF_runs','Impact Location by Time', flat_lpfs)
#print(noiselist)
makegraphpg(scriptsdir,'noise_runs', 'Noise Location by Time',noise)
makegraphpg(scriptsdir,'momenta_graphs_runs', 'Momenta Distribution by Time',momenta_graphs)
makegraphpg(scriptsdir,'skyloc_runs', 'Skylocation by Time',skyloc)
makegraphpg(scriptsdir,'likelihood_runs','Likelihood Chain by Time',logprob)
makegraphpg(scriptsdir,'impact_runs','Impact Model by Time',impacts)
#makegraphpg(homedir,'error_dofs_runs', dofs)

mom_filenames = np.asarray(mom_filenames)
#makelistpage(webaddress,filenames, title, capital, scriptsdir):
makelistpage(webaddress = webaddress,filenames = mom_filenames, title = 'momenta',search ='', capital = 'Momenta', scriptsdir = scriptsdir)
makelistpage(webaddress = webaddress,filenames = run_filenames,title ='time', search = nameruns, capital = 'Time',scriptsdir = scriptsdir)
