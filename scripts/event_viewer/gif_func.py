import imageio
import re
import glob
import os
from os import listdir
from os.path import isfile, join

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)



def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

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
#homedir = '/home/sophie.hourihane/REU'
##'/home/sophie.hourihane/public_html/website/public_html/gifs'
#public_runs = homedir + '/runs' #'/home/sophie.hourihane/public_html/website/public_html/images/runs/'
#
##public_moms = '/home/sophie.hourihane/public_html/website/public_html/images/momenta/'
#
def gif_maker(gifdir, currun, mednmom, savedir):
#	mom_rundirs,gifdirnames = get_the_subdir(gifdir)
#	for mom_run in mom_rundirs:
#		#getnames = mom_run.replace(gifdir,"")
#		#currun = getnames[-10:]
#		#mednmom = getnames[5:17]
#		#print(currun, mednmom)
	topbotdirs,namemomrun = get_the_subdir(gifdir)
	for topbot in topbotdirs:
		view = topbot[-3:]
		print(view)
		filenames = [fil for fil in listdir('%s'%(topbot)) if isfile(join('%s'%(topbot),fil))]
		images = []
		filenames.sort(key=alphanum_key)
		for filename in filenames:	
			images.append(imageio.imread(topbot+'/'+filename))
		
		imageio.mimsave('%s/gif_run%s_mom%s_%s.gif'%(savedir,currun,mednmom,view), images)
		
		#imageio.mimsave('%s/mom_%s/gif_run%s_mom%s_%s.gif'%(public_moms,mednmom,currun,mednmom,view), images)
