'''code written by Matt Wilde to efficiently drizzle already tweaked images suitable for 
	photometry and blinking. Must run in the /Users/mattwilde/47Tuc/BLINK/ ... with 
	appropriate detectrfilter
'''



from glob import glob
#from astropy.io import fits
from drizzlepac import *
from stsci.tools import teal
teal.unlearn('updatenpol')
teal.unlearn('astrodrizzle')


'''
need to implement this to be extensible to other directories:
updatenpol.update(list_of_images,refdir='/path/to/reference/files/',\
          local=True)
'''


###seams to only work in python, not pyraf...
def drizzler(filter,detector):
	'''function to auto generate drizzle functions based on a given filter and detector. 
	Takes in a filter and detector as arguments in quotations seperated by a comma.
	Ex: drizzler('F435W','WFPC2')
	'''
	
	epochs1 = glob('9*')
	epochs2 = glob('0*')
	epochs3 = glob('12*')
	epochs4 = glob('13*')
	epochs = epochs1 + epochs2 + epochs3 + epochs4
	
	
	if detector == 'WFPC2':
	
		for e in epochs:
			if 'single' in e:
				print 'single input image: go back and do separately' 
			else:
				astrodrizzle.AstroDrizzle(input='@'+e, driz_sep_fillval='99999', \
					output='O_PHOT_PC_'+detector+'_'+filter+'_'+e, wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
					updatewcs='False', driz_sep_bits='8,1024',combine_type='minmed', \
					driz_cr_snr='10 8.5', driz_cr_scale='6.0 5.5', final_wht_type='EXP', \
					final_wt_scl='exptime', final_pixfrac='1.0', final_fillval='99999', \
					driz_cr_corr= 'True', final_bits='8,1024',group='1',final_units='counts')
		
				astrodrizzle.AstroDrizzle(input='@'+e, driz_sep_fillval='99999', \
					output='O_PHOT_WF2_'+detector+'_'+filter+'_'+e, wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
					updatewcs='False', driz_sep_bits='8,1024',combine_type='minmed', \
					driz_cr_snr='7.0 5.5', driz_cr_scale='2.0 1.5', final_wht_type='EXP', \
					final_wt_scl='exptime', final_fillval='99999', \
					driz_sep_wcs='True', driz_sep_scale=0.1, final_wcs='True',final_scale=0.1, \
					driz_cr_corr= 'True', final_bits='8,1024',group='2',final_units='counts')  
		
				astrodrizzle.AstroDrizzle(input='@'+e, driz_sep_fillval='99999', \
					output='O_PHOT_WF3_'+detector+'_'+filter+'_'+e, wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
					updatewcs='False', driz_sep_bits='8,1024',combine_type='minmed', \
					driz_cr_snr='7.0 5.5', driz_cr_scale='2.0 1.5', final_wht_type='EXP', \
					final_wt_scl='exptime', final_fillval='99999', \
					driz_sep_wcs='True', driz_sep_scale=0.1, final_wcs='True',final_scale=0.1, \
					driz_cr_corr= 'True', final_bits='8,1024',group='3',final_units='counts')    
		
				astrodrizzle.AstroDrizzle(input='@'+e, driz_sep_fillval='99999', \
					output='O_PHOT_WF4_'+detector+'_'+filter+'_'+e, wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
					updatewcs='False', driz_sep_bits='8,1024',combine_type='minmed', \
					driz_cr_snr='7.0 5.5', driz_cr_scale='2.0 1.5', final_wht_type='EXP', \
					final_wt_scl='exptime', final_fillval='99999', \
					driz_sep_wcs='True', driz_sep_scale=0.1, final_wcs='True',final_scale=0.1, \
					driz_cr_corr='True', final_bits='8,1024',group='4',final_units='counts')
		
				astrodrizzle.AstroDrizzle('@'+e, driz_sep_fillval='99999', \
					output='O_BLINK_'+detector+'_'+filter+'_'+e, wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
					updatewcs='False', driz_sep_bits='8,1024',combine_type='minmed', \
					driz_cr_snr='2.5 2.0', driz_cr_scale='1.0 0.6', final_wht_type='EXP', \
					final_wt_scl='exptime', final_pixfrac='1.0', final_fillval='99999', \
					driz_cr_corr='True', final_bits='8,1024',final_units='cps')
		
		
	if detector == 'WFC':
		for e in epochs:
			if 'single' in e:
				print 'single input image: go back and do separately' 
				break
			astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
				output='O_BLINK_'+detector+'_'+filter+'_'+e, \
				combine_type='minmed',driz_cr='True',final_pixfrac='1.0',final_units='cps')
		
			astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+detector+'_'+filter+'_PC', \
				output='O_PHOT_'+detector+'_'+filter+'_'+e, \
				combine_type='minmed',driz_cr='True',final_pixfrac='1.0',final_units='counts')
		
			#print 'TWEAK_'+detector+'_'+filter+'_PC'

	if detector == 'WFC3':
		for e in epochs:
			if 'single' in e:
				print 'single input image: go back and do separately' 
				break
			astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+detector+'_'+filter+'_PC_1', \
				output='O_BLINK_'+detector+'_'+filter+'_'+e, \
				combine_type=' minmed',driz_cr='True',final_pixfrac='1.0',final_units='cps')
		
			astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+detector+'_'+filter+'_PC_1', \
				output='O_PHOT_'+detector+'_'+filter+'_'+e, \
				combine_type='minmed',driz_cr='True',final_pixfrac='1.0',final_units='counts')
						
	if detector == 'HRC':
		#epochs=['050302_single']
		for e in epochs:
			if 'single' in e:
				e = e[:-7]
				astrodrizzle.AstroDrizzle(input='@'+e+'_single', \
				output='O_BLINK_'+detector+'_'+filter+'_'+e,final_units='cps', \
				configobj='single.cfg')
				
				astrodrizzle.AstroDrizzle(input='@'+e+'_single', \
				output='O_PHOT_'+detector+'_'+filter+'_'+e,final_units='counts', \
				configobj='single.cfg')
				
				
 
			else:
				astrodrizzle.AstroDrizzle(input='@'+e, \
					output='O_BLINK_'+detector+'_'+filter+'_'+e,final_units='cps', \
					configobj='multi_cps.cfg')
				
				astrodrizzle.AstroDrizzle(input='@'+e, \
					output='O_PHOT_'+detector+'_'+filter+'_'+e,final_units='counts', \
					configobj='multi_counts.cfg')
'''				
				astrodrizzle.AstroDrizzle(input='@'+e, \
					output='O_BLINK_'+detector+'_'+filter+'_'+e, \
					combine_type=' minmed',driz_cr='True',final_pixfrac='1.0',final_units='cps')
	
				astrodrizzle.AstroDrizzle(input='@'+e, \
					output='O_PHOT_'+detector+'_'+filter+'_'+e, \
					combine_type='minmed',driz_cr='True',final_pixfrac='1.0',final_units='counts')
'''


'''		
		else:
			for e in epochs:
				if 'single' in e:
					print 'single input image: go back and do separately' 
					break
				astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+filter+'_PC', \
					output='O_BLINK_'+detector+'_'+filter+'_'+e, \
					combine_type=' minmed',driz_cr='True',final_pixfrac='1.0',final_units='cps')
		
				astrodrizzle.AstroDrizzle(input='@'+e,wcskey='TWEAK_'+filter+'_PC', \
					output='O_PHOT_'+detector+'_'+filter+'_'+e, \
					combine_type='minmed',driz_cr='True',final_pixfrac='1.0',final_units='counts')
'''