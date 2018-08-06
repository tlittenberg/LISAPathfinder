import os
import re
import numpy as np
import microTools as mT
from impactClass import impactClass
import pathlib

##### Constants ####
BASE_DIR = "/Users/shouriha/LISAPathfinder/scripts/Analysis" ##"/home/sophie.hourihane/public_html/viterbi"
WEB_ADDRESS = BASE_DIR #"https://ldas-jobs.ligo-wa.caltech.edu/~sophie.hourihane/viterbi" 
img_directory = BASE_DIR + "/data"

################# COLORS ###############
#complementary
#COLORS = ["#D04B3F", "#487384", "#CC8F33", "#2C2112"]
#green
COLORS = ["#EBF7E3","#9BD770","#66B032","#375F1B","#1B3409"]

################# STYLE #################
START = '''
<!DOCTYPE html>
<html>
<link rel="stylesheet" href="https://www.w3schools.com/w3css/4/w3.css">
<head>
'''

STYLE = '''
<style>
body {
	background-color:%s; '''%COLORS[0] + '''
	float: center;
}
h2 {
	text-align:center;
}
div.gallery {
	margin: 5px;
	border: 1px solid #ccc;
	float: left;
	width: 600px;
}
div.gallery:hover {
	border: 1px solid #777;
}
div.gallery img {
	width: 100%%;
	height: auto;
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
	background-color: %s;''' %(COLORS[0]) + '''

	color: black;
}

.next {
	background-color: %s;''' %(COLORS[len(COLORS) - 1]) + '''
	color: white;
}

.round {
	border-radius: 50%%;
}



table, th, td{
	border: 1px #fff;
	padding: 15px;
	text-align:left;
	margin-left: auto;
	margin-right: auto;
}  ''' + '''
th {
	background-color:%s; '''%(COLORS[len(COLORS) - 1]) + '''
	color:%s;'''%(COLORS[0]) + '''
}
th:hover{
	background-color:#fff;
	color:%s;'''%(COLORS[len(COLORS) - 1]) + '''
}
tr:nth-child(even) { '''+ '''
	background-color: %s; '''%(COLORS[0]) + '''
}
tr:hover{
	background-color:#fff;
}
table#params{
	width: 70%%;''' + '''
	background-color: %s;'''%(COLORS[1]) + '''
}
table#header{
	float: center;
	width: 70%%;''' + '''
	background-color: %s;'''%(COLORS[0]) + '''
}


ul {
	list-style-type: none;
	margin: 0;
	padding: 0;
	overflow: hidden; ''' + '''
	background-color: %s;'''%(COLORS[len(COLORS) - 1]) + '''
}

li {
	float: left;
}

li a, .dropbtn {
	display: inline-block;
	color: white;
	text-align: center;
	padding: 14px 16px;
	text-decoration: none;
}

li a:hover, .dropdown:hover .dropbtn { ''' + '''


	background-color: %s'''%(str(COLORS[3])) + ''';
}


li.dropdown {
	display: inline-block;
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
	color: black;
	padding: 12px 16px;
	text-decoration: none;
	display: block;
	text-align: left;
}

.dropdown-content a:hover {background-color: #f1f1f1}

.dropdown:hover .dropdown-content {
	display: block;
}
</style> '''

HEADER = '''

<ul>
	<li><a href="index.html">Home</a></li>
	<li class="dropdown">
	  <a href="javascript:void(0)" class="dropbtn">More Plots</a>
	  <div class="dropdown-content">
		<a href="level_curves.html">Level Curves</a>
		<a href="#">Link 2</a>
		<a href="#">Link 3</a>
	  </div>
	</li>
</ul>


</style> '''

FUNCTIONS = '''
<script>

function sortTable(n) {
	var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
	table = document.getElementById("params");
	switching = true;
	//Set the sorting direction to ascending:
	dir = "asc"; 
	/*Make a loop that will continue until
	no switching has been done:*/
	while (switching) {
	  //start by saying: no switching is done:
	  switching = false;
	  rows = table.getElementsByTagName("TR");
	  /*Loop through all table rows (except the
	  first, which contains table headers):*/
	  for (i = 1; i < (rows.length - 1); i++) {
		//start by saying there should be no switching:
		shouldSwitch = false;
		/*Get the two elements you want to compare,
		one from current row and one from the next:*/
		x = rows[i].getElementsByTagName("TD")[n];
		y = rows[i + 1].getElementsByTagName("TD")[n];
		/*check if the two rows should switch place,
		based on the direction, asc or desc:*/
		if (dir == "asc") {
		  if (Number(x.innerHTML.toLowerCase()) > Number(y.innerHTML.toLowerCase())) {
			//if so, mark as a switch and break the loop:
			shouldSwitch= true;
			break;
		  }
		} else if (dir == "desc") {
		  if (Number(x.innerHTML.toLowerCase()) < Number(y.innerHTML.toLowerCase())) {
			//if so, mark as a switch and break the loop:
			shouldSwitch = true;
			break;
		  }
		}
	  }
	  if (shouldSwitch) {
		/*If a switch has been marked, make the switch
		and mark that a switch has been done:*/
		rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
		switching = true;
		//Each time a switch is done, increase this count by 1:
		switchcount ++;      
	  } else {
		/*If no switching has been done AND the direction is "asc",
		set the direction to "desc" and run the while loop again.*/
		if (switchcount == 0 && dir == "asc") {
		  dir = "desc";
		  switching = true;
		}
	  }
	}
}



</script>'''



def getFiles(directory):
	files = []
	try:
		for fil in os.listdir(directory):
			if os.path.isfile(os.path.join((directory),fil)):
				files.append(os.path.join((directory),fil))
		return files
	except FileNotFoundError:
		print("Not Found", directory)
		return []

def getDirs(in_path):
	subdirs = []
	for x in os.walk(in_path):
		subdirs.append(x[0])
	return subdirs

#takes in absolute path and returns webaddress
def makeSRC(filename):
	left_over = filename.replace(BASE_DIR, "")
	src = WEB_ADDRESS + left_over
	return src

#takes in list of (absolute) files, creates an image gallery
def makeGallery(onlyfiles, name = 'TEST'):
	gallery_string = ""
	for i, f in enumerate(onlyfiles):
		gallery_string += '''<div class = "gallery %s">\n'''%(name)
		#create image gallery
		img_src = makeSRC(f)
		img_formatted = '''\t <img src ="''' + img_src + '''" alt ="''' + f +'''">\n''' 
		gallery_string += img_formatted
		gallery_string += "</div>\n"
	return gallery_string

def makeParamTable(params, ignore = "VOID", class_name = "params", link = "VOID"):
	# Makes params list of params
	if type(params) is not list: params = [params]
	if type(ignore) is not list: ignore = [ignore]
	
	#had vary table id or else sort doesnt work
	table_script = '''<table id="%s"> \n'''%(class_name)
	all_cols = ['isValid', 'segment', 'gps', 'N', 'snr', 't0', 'Ptot', 
				'lat', 'lon', 'rx', 'ry', 'rz', 'face']

	headers = []
	data = []
	# Make header
	count = 0
	for i, col in enumerate(all_cols):
		if i == 0:
			table_script += "\t<tr>\n"
		if col in ignore:
			continue
		else:
			table_script += '''\t\t<th onclick="sortTable(%i)">'''%(count) + col + "</th> \n"
			count += 1
		if i == len(all_cols) - 1:
			table_script += "\t</tr>\n"
	# make Data
	for j, param in enumerate(params):
		for i, col in enumerate(all_cols):
			if i == 0:
				# links to webpage
				if link == 'page':
					table_script += '''\t<tr onclick="window.location='%s';">\n'''%((WEB_ADDRESS + 
													"/html_scripts/" + param.filename() + ".html")) 
				else:
					table_script += '''\t<tr>'''
			if col in ignore:
				continue
			else:
				if col == 'Ptot' or col == 't0':
					table_script += "\t\t<td>" + "%1.3e </td> \n" % param.getMedian(col) 
				else:
					table_script += "\t\t<td>" + "%.3f </td> \n" % param.getMedian(col) 
			if i == len(all_cols) - 1:
				table_script += "\t</tr>\n"
	#Close table
	table_script += "</table>\n"
	return table_script
			
def writeScript(filename, body):
	script = START + STYLE + "\n</head>"
	script += "\n<body>"
	script += HEADER
	script += body
	script += FUNCTIONS
	script += "\n</body>"
	script += "\n</html>"

	f = open(BASE_DIR + "/html_scripts/" + filename, "w+")
	f.write(script)
	f.close()

#returns correctly named images from given parameter
def findParamPlots(param, find = 'VOID', ignore = "VOID", directory = BASE_DIR + "/plots/"):
	if type(ignore) is not list: ignore = [ignore]
	if type(find) is not list: find = [find]
	images = []
	only_files = getFiles(directory)
	for f in only_files:
		# only images
	
		if not ".png" and not '.gif' in f:
			continue

		# If we are looking for specific files
		if find != 'VOID':
			for fin in find:
				if fin in f:
					images.append(f)
				else:
					continue
		else:
			images.append(f)
	return images
	

def makeParamPage(params, next_param, prev_param):
	print(params.filename(), "Make Param Page")
	print(BASE_DIR)
	filename = params.filename() + ".html"
	table = makeParamTable(params, class_name = "header", ignore = ["face"], link = 'page' )

	p = pathlib.PurePath(os.getcwd())
	baseDir = str(p.parent)
	dataPath = pathlib.Path(baseDir + '/data')
	pickles = list(dataPath.glob('*_grs1.pickle'))

	galleries = ''

	#find images
	image_list = findParamPlots(params, find = ['dual', 'HTC_JFC_pop'], ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list, name = 'dual') + '<br>'
	
	image_list = findParamPlots(params, find = ['pop'], ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list, name = 'dual') + '<br>'

	image_list = findParamPlots(params, find = 'sky', ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list)

	image_list = findParamPlots(params, find = 'sun', ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list)

	image_list = findParamPlots(params, find = 'flat', ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list)

	image_list = findParamPlots(params, find = '.gif', ignore = ['face'], directory = BASE_DIR + '/plots/' + params.filename())
	galleries  += makeGallery(image_list)

	# Makes Body of Page 
	body = """ <a href="%s" class="previous">&laquo; Previous</a> \n"""%(prev_param.filename() + '.html')
	body += """<a href="%s" class="next">Next &raquo;</a> \n"""%(next_param.filename() + '.html')
	body += "<h2> Impact Parameters</h2>\n"
	body += table
	body += "<h2> Impact Images </h2>"
	body += galleries

	writeScript(filename, body)

def makeLevelCurves():
	#get only level curve files
	level_curves = []
	h_val = []
	only_files = getFiles(BASE_DIR + "/plots")
	for f in only_files:
		if "level_curves" in f:
			level_curves.append(f)
		if "h_val" in f:
			h_val.append(f)
	gallery = makeGallery(level_curves)
	gallery_h_val = makeGallery(h_val)

	body = gallery + gallery_h_val

	writeScript("level_curves.html", body)
	

def makeIndex():

	p = pathlib.PurePath(os.getcwd())
	baseDir = str(p.parent)
	dataPath = pathlib.Path(baseDir + '/data')
	pickles = list(dataPath.glob('*_grs1.pickle'))

	onlyfiles = getFiles(BASE_DIR + "/plots")

	#Lists every family of max_viterbi runs i've done
	directories = getDirs(BASE_DIR + '/data')
	param_list = []


	# run through the pickles
	segments = []
	for p in pickles:
		segments.append(str(p.stem[0:10]))
	segments = np.sort(segments)

	i = -1
	for s in segments:
		i += 1
		# identify segment
		segment = s#str(p.stem[0:10])


		# load GRS1 data
		chainFile = baseDir + '/data/' + str(segment) +'_grs1.pickle'
		try:
			param = impactClass(chainFile)
			if not os.path.isdir(baseDir + '/plots/' + param.filename()):
				continue

		except ValueError:
			continue
		else:
			# Make Links to Next GPS Time
			# Get previous segment
			if i == 0:
				prev_seg = segment
			else:
				prev_seg = segments[i - 1]

			if i == len(pickles) - 1:
				next_seg = segment
			else:
				next_seg = segments[i + 1]

			chainFilenext = baseDir + '/data/' + str(next_seg) +'_grs1.pickle'
			chainFileprev = baseDir + '/data/' + str(prev_seg) +'_grs1.pickle'
			next_param = impactClass(chainFilenext)
			prev_param = impactClass(chainFileprev)

			makeParamPage(param, next_param, prev_param)
			param_list.append(param)
			

	#Links to Impact pages
	table = makeParamTable(param_list, ["face"], link =  'page')
	print(table)

	body = ""
	body += '''<h2>List of Impacts</h2> \n''' + table

	writeScript("index.html", body)

makeIndex()
