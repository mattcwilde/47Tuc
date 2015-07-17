'''This program reads in the photometry information produced by daophot on the 
	CVs V1, V2, AKO9, V3 in 47 Tucanae. It will also read in the epoch/time 
	information from each images header and plot magnitude v epoch. 
	
	!!!If i were to do this again I would make separate programs, the first being one
	that made a table out of the data (also creating lateX file)
	
	The other programs would plot the data'''

from math import *
#import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob
from astropy.io import fits
import pylab as P

P.close('all')

tab = Table()

def NewMag(mag,pflam,expt):
	"take measured magnitudes and convert to STMAGS using zpt"

	z = -2.5*log10(pflam) - 21.1
	new = mag - 25 + z + 2.5*log10(expt)
	result = new
	return result

def limits(object):
	if object == 'V1':
		P.ylim(24,18)
	elif object == 'V2':
		P.ylim(23,18)
	elif object == 'AKO9':
		P.ylim(21,16)
	elif object == 'V3':
		P.ylim(30,13)
	else:
		print "wrong object"

def filter_plot(data,merr,time,filter_name,object_name,bad,bad_t):
	"create plots of individaual fiters"
	fig = P.figure()
	ax = fig.add_axes([0.1,0.1,0.8,0.8])
	P.scatter(time,data,color='k')
	P.scatter(bad_t,bad,marker=None)
	if len(bad_t) != 0:
		ax.errorbar(bad_t, bad, yerr=0, lolims=True,capsize=6,linestyle='None',color='k')
	if len(data) != 0:
		ax.errorbar(time, data, yerr=merr,linestyle='None', color='k')
	ax.set_xlabel('MJD')
	ax.set_ylabel('{0} {1} (STMAGNITUDE)'.format(object_name, filter_name))
	fig.gca().invert_yaxis()
	P.xlim(48960,56400)
	#P.ylim(30,13)
	limits(object_name) #create ylims using function
	fig.savefig('{0}_{1}_timeseries.png'.format(object_name,filter_name))
	#P.show()
	P.close()
	


def color_plot(d,object):
    "make plot of all filters on one plot, differentiated by color"
    
    fig = P.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    P.title(object)
    ax.set_xlabel('MJD')
    ax.set_ylabel('{} (STMAGNITUDE)'.format(object))
    limits(object)
    #P.ylim(,13)
    P.xlim(48960,59400)
    cols = ['#0000FF','#4000FF','#0040FF','#0080FF','#00BFFF','#00FFFF','#00FFBF',
    	'#00FF80','#00FF40','#40FF00','#BFFF00','#FFFF00','#FFBF00',
    	'#FF8000','#FF4000','#FF0000']
    for k,c in zip(sorted(d.keys()),cols):
	    t = d[k][2]
	    y = d[k][0]
	    bad_t = d[k][4]
	    bad = d[k][3]
	    err = d[k][1]
	    P.plot(t,y,'o',c=c,label=k)
	    P.scatter(bad_t,bad,marker=None)
	    if len(bad)!= 0:
	    	print object
	    	ax.errorbar(d[k][4],d[k][3], yerr=0, lolims=True, capsize=6,linestyle='None')
	    	#ax.errorbar(d[k][4],d[k][3], yerr=0, lolims=True, capsize=6)
	    else:
	    	ax.errorbar(t, y, yerr=err,linestyle='None', color='black')
	    #ax.errorbar(t,y,yerr=err,linestyle='None',c=c)
	    
	    P.legend()
	#return None
	
	
#initialize arrays to store each CVs measured data
objects = ['V1','V2','AKO9','V3']

vio = ['F225W', 'F250W', 'F300W', 'F330W', 'F336W']
blue = ['F410M', 'F435W', 'F439W', 'F450W', 'F475W']
vis = ['F550M', 'F555W', 'F606W']
red = ['F625W', 'F775W', 'F814W']

filtd_v1 = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}

filtd_v2 = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}

filtd_ako9 = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}
	
filtd_v1_bad = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}

filtd_v2_bad = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}

filtd_ako9_bad = {'F225W': [],'F250W': [], 'F300W': [], 'F330W': [], 'F336W': [], 
	'F410M': [], 'F435W': [], 'F439W': [], 'F450W': [], 'F475W': [],
	'F550M': [], 'F555W': [], 'F606W': [], 'F625W': [], 'F775W': [], 'F814W': []}
	
filters = vio+blue+vis+red
#filters = ['F439W']
''' for each mag file read in mag and merr for each object'''

for f in filters:
	P.clf()
	file_list = glob('O*cnts*'+f+'*mag.1')
	#v3_file_list = glob('*mag.2')
	#file_list = glob('O_small_cnts_WFPC2_F439W_950901_drz_sci.fits.mag.1')

	v1_phot = []
	v1_phot_err = []
	v1_shift = []
	v1_bad = []
	v1_bad_time = []
	v1_bad_err = []

	v2_phot = []
	v2_phot_err = []
	v2_shift = []
	v2_bad = []
	v2_bad_time = []
	v2_bad_err = []

	ako9_phot = []
	ako9_phot_err = []
	ako9_shift = []
	ako9_bad = []
	ako9_bad_time = []
	ako9_bad_err = []

	filt_list = []
	name_list = []

	epoch = []

	v1_epoch = []
	v2_epoch = []
	ako9_epoch = []

	#Make functions to correct mags to STMAG

	#read through mag files and fits files to get the objects magnitudes,
	#magnitude errors, as well as photflam and exposure time. 
	for i,file in enumerate(file_list):
	
		#read daophot file
		t = Table.read(file,format='ascii.daophot')
	
		#find root of file name (remove .mag.1)
		im_name = file[:-6]
	
		#read in epoch information (JD)
		hdulist = fits.open(im_name)
		photflam = hdulist[0].header['PHOTFLAM']
		expt = hdulist[0].header['EXPTIME']
		tau = hdulist[0].header['EXPSTART']
		epoch.append(tau)
	
		v1_epoch.append(tau)
		v2_epoch.append(tau)
		ako9_epoch.append(tau)

		#get filter name from image name (get detector too??)
		filt_index = im_name.find('_F')
		filt = (im_name[filt_index+1:filt_index+6])
		filt_list.append(filt)
		name_list.append(im_name)
	
		v1_shift = sqrt(t['XSHIFT'][0]**2 + t['YSHIFT'][0]**2)
		v2_shift = sqrt(t['XSHIFT'][1]**2 + t['YSHIFT'][1]**2)
		ako9_shift = sqrt(t['XSHIFT'][2]**2 + t['YSHIFT'][2]**2)
	
		v1_error_code = t['CERROR'][0]
		v2_error_code = t['CERROR'][1]
		ako9_error_code = t['CERROR'][2]
		
		if v1_error_code == 'LowSnr':
			print 'v1 not on image:',im_name
		if v2_error_code == 'LowSnr':
			print 'v2 not on image:',im_name
		if ako9_error_code == 'LowSnr':
			print 'ako9 not on image:',im_name
		
		
		#fill the arrays with phot data applying STMAG correction
		if v1_shift > 4:
			v1_bad.append(NewMag(22.013,photflam,expt))
			v1_bad_err.append(0)
			t1 = v1_epoch.pop() #remove the bad time from the list of good ones
			v1_bad_time.append(t1) #put bad time into list 
		else:
			v1_phot.append(NewMag(t['MAG'][0],photflam,expt))
			v1_phot_err.append(t['MERR'][0])
		
		if v2_shift > 4:
			v2_bad.append(NewMag(22.013,photflam,expt))
			v2_bad_err.append(0)
			t1 = v2_epoch.pop()
			v2_bad_time.append(t1)
		else:
			v2_phot.append(NewMag(t['MAG'][1],photflam,expt))
			v2_phot_err.append(t['MERR'][1])
		
		if ako9_shift > 4:
			ako9_bad.append(NewMag(22.013,photflam,expt))
			ako9_bad_err.append(0)
			t1 = ako9_epoch.pop()
			ako9_bad_time.append(t1)
		else:
			ako9_phot.append(NewMag(t['MAG'][2],photflam,expt))
			ako9_phot_err.append(t['MERR'][2])	
		
		hdulist.close()
	'''		
	set_filt_list = list(set(filt_list))
	set_filt_list.sort()
	if sorted(set_filt_list) != filters:
		print "set_filt_list != filters"
	'''
	
	filter_plot(v1_phot,v1_phot_err,v1_epoch,f,'V1',v1_bad,v1_bad_time)
	filter_plot(v2_phot,v2_phot_err,v2_epoch,f,'V2',v2_bad,v2_bad_time)
	filter_plot(ako9_phot,ako9_phot_err,ako9_epoch, f, 'AKO9',ako9_bad,ako9_bad_time)
	
	###make dictionary containing photometry, error, and epoc
	filtd_v1[f] = [v1_phot,v1_phot_err,v1_epoch,v1_bad,v1_bad_time]
	filtd_v2[f] = [v2_phot,v2_phot_err,v2_epoch,v2_bad,v2_bad_time]
	filtd_ako9[f] = [ako9_phot,ako9_phot_err,ako9_epoch,ako9_bad,ako9_bad_time]
	
	#diagnostic
	if len(v2_bad) != 0:
		print 'object too faint in',im_name
	
	
color_plot(filtd_v1,'V1')
color_plot(filtd_v2,'V2')
color_plot(filtd_ako9,'AKO9')
