from astropy.io import fits
import glob,os,shutil
import sys,csv
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.ndimage import zoom
#from numpy import *
#from gaussian_model_img import *
from casatools import table
import random
import bdsf
from scipy.optimize import curve_fit
import pylab as pl
from matplotlib.ticker import AutoLocator
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import *
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from casatasks.private.simutil import simutil
import casatools
from casatools import simulator
from casatools import measures
from casatasks import importfits, tclean, exportfits, simobserve,gaincal,applycal,mstransform
from casatools import simulator, measures, quanta, table
from casaplotms import plotms

sm = simulator()
me = measures()
qa = quanta()
tb = table()

def header(image,RA,DEC,delta):
 imagename=image[:-5]+'_'+str(RA)+'deg'+'_'+str(DEC)+'deg'+'.fits'
 hdu=fits.open(image)[0]
 data=hdu.data
 data=data#*1e6
 header=hdu.header
 if len(data.shape)==4:
  (a,b,xx,yy)=data.shape
 if len(data.shape)==3:
  (a,xx,yy)=data.shape
 if len(data.shape)==2:
  (xx,yy)=data.shape
 xc=xx/2.;yc=yy/2.
 RA=RA
 DEC=DEC
 cdelt1=-(delta/3600.)
 cdelt2=(delta/3600.)      
 header['BSCALE'] =  1.                                                  
 header['BZERO'] =   0.                                                  
 header['BUNIT'] = 'Jy/pixel'                           
 header['BMAJ'] =  0.00333886280458654                                                  
 header['BMIN'] =  0.00304123673936456                                                 
 header['BPA'] =   129.347288653168                                                
 header['EQUINOX'] = 2000.                                      
 header['BTYPE'] = 'Intensity'                                                         
 header['CTYPE1'] = 'RA---SIN'                          
 header['CRPIX1'] =  xc                                                  
 header['CRVAL1'] = RA         # RA                                         
 header['CDELT1'] = cdelt1                                                
 header['CUNIT1'] = 'deg     '                                                            
 header['CTYPE2'] = 'DEC--SIN'                              
 header['CRPIX2'] =  yc                                                  
 header['CRVAL2'] =  DEC                 #Dec                                    
 header['CDELT2'] =  cdelt2                                                 
 header['CUNIT2'] = 'deg     '                                                            
 header['CTYPE3'] = 'FREQ    '                                        
 header['CRPIX3'] = 1.                                                  
 header['CRVAL3'] = 2500e6                                               
 header['CDELT3'] = 855791015.625                                                  
 header['CUNIT3'] = 'Hz      '                                                            
 header['CTYPE4'] = 'STOKES  '                                                            
 header['CRPIX4'] = 1.                                                  
 header['CRVAL4'] = 1.                                                  
 header['CDELT4'] = 1.                                                  
 header['CUNIT4'] ='        '                                                 
 
 img=fits.writeto(imagename,data,header,overwrite=True)
 return img


def rescale_galactic_center_image_to_andromeda(
    input_fits,
    output_fits,
    gc_distance_kpc=8.2,
    m31_distance_kpc=780,
    scale_brightness=False
):
    """
    Rescale a Galactic Center FITS image to simulate how it would appear at the distance of Andromeda (M31).
    
    Parameters:
    ----------
    input_fits : str
        Path to the input FITS file (MeerKAT GC image).
    output_fits : str
        Path to save the rescaled output FITS file.
    gc_distance_kpc : float
        Distance to the Galactic Center (default = 8.2 kpc).
    m31_distance_kpc : float
        Distance to M31 (default = 780 kpc).
    scale_brightness : bool
        Whether to scale brightness by inverse square of distance ratio.
        Only use True if input is in Jy/pixel or total flux.
    """
    # Load image
    hdul = fits.open(input_fits)
    data = hdul[0].data.squeeze()
    header = hdul[0].header
    wcs = WCS(header)

    # Compute scale
    scale_factor = gc_distance_kpc / m31_distance_kpc
    flux_factor = scale_factor ** 2
    print(f"[INFO] Angular scale factor: {scale_factor:.4f}")
    print(f"[INFO] Flux scale factor (if applied): {flux_factor:.2e}")

    # Resize image
    zoomed_data = zoom(data, scale_factor, order=3)  # cubic interpolation

    # Brightness scaling (optional)
    if scale_brightness:
        zoomed_data *= flux_factor

    # Adjust header
    header['CDELT1'] = header['CDELT1'] / scale_factor
    header['CDELT2'] = header['CDELT2'] / scale_factor

    ny, nx = zoomed_data.shape
    header['NAXIS1'] = nx
    header['NAXIS2'] = ny
    header['CRPIX1'] = nx / 2
    header['CRPIX2'] = ny / 2
    header['CRVAL1'] = 10.68470833  # RA of M31 in deg
    header['CRVAL2'] = 41.269065    # Dec of M31 in deg
    header['OBJECT'] = 'GC_model_at_M31'

    # Save new FITS
    hdu = fits.PrimaryHDU(zoomed_data, header=header)
    hdu.writeto(output_fits, overwrite=True)
    print(f"[SUCCESS] Saved rescaled image to: {output_fits}")



def psf_size(vis):
  ms.open(vis)
  ms.selectinit(datadescid=0)
  ms.selectinit(reset=True)
  uvd=ms.range(["uvdist"])
  uvdistance=uvd['uvdist'][1]
  spwInfo = ms.getspectralwindowinfo()
  tb.open(vis + "/SPECTRAL_WINDOW")
  chan_freq = tb.getcol("CHAN_FREQ")
  tb.close()
  tb.done()
  distance=uvdistance
  mean_freq=mean(chan_freq)
  c=3e8 #in cm ;3e8 in m
  psf_size=((c/mean_freq)/distance) *(180/3.14)*3600
  psf_size1=psf_size/3.;psf_size2=psf_size/4.
  return (distance,psf_size1,psf_size2)

def fov(vis):
  ms.open(vis)
  ms.selectinit(datadescid=0)
  ms.selectinit(reset=True)
  spwInfo = ms.getspectralwindowinfo()
  tb.open(vis + "/SPECTRAL_WINDOW")
  chan_freq = tb.getcol("CHAN_FREQ")
  tb.close()
  tb.done()
  mean_freq=mean(chan_freq)
  PB=42.0 * 1e9/mean_freq
  return (PB)

def generateOffsetPerAntennaPerTime(msname: str ="", maxoffset: str ="1arcmin", distribution: str ="uniform", sigma: str ="1arcmin", ):
	"""
	distribution can be uniform or gaussian .gaussian need sigma, uniform need maxoffset
	return dictionary with key (ant, time) and pair of x,y offset
	"""
	tb=casatools.table()
	tb.open(f"{msname}/ANTENNA")
	nant=tb.nrows()
	tb.close()
	antids=np.arange(nant)
	tb.open(msname)
	t=np.unique(tb.getcol("TIME"))
	tb.done()
	maxoff=qa.convert(qa.quantity(maxoffset), 'rad')['value']
	sig=qa.convert(qa.quantity(sigma), 'rad')['value']
	offsetdict={}
	for tim in t:
		for ant in antids:
			key=(ant, tim)
			ang=random.uniform(0, 2*np.pi)
			ampoff=random.uniform(0.0, maxoff)
			offset=[ampoff*np.cos(ang), ampoff*np.sin(ang)]
			offsetdict[key]=offset
	return offsetdict
	
def distFromCent(msname: str ="", field : int =0, compdirection={}):
	"""
	returns a pair of distance from phasecenter in radians
	"""
	msmd.open(msname)
	ph=msmd.phasecenter(field)
	msmd.done()
	print(f"{ph} , {compdirection}, {me.separation(ph, compdirection)}")
	dist=qa.convert(me.separation(ph, compdirection), "rad")['value']
	pa=qa.convert(me.posangle(ph, compdirection), 'rad')['value']
	val=[dist*np.sin(pa), dist*np.cos(pa)]
	return val


def simulatePontingError(mode: str = "single", msname: str = "", nchan: int = 1,
						 centrefreq: str = "93GHz",chanwidth: str = "10MHz",
						 inttime: str = "10s", obstime: str = "1s",  noise: str = "0.1mJy",
						 numpointsource: int =10, maxflux: str = "1Jy", maxpointError: str = "5arcsec"):
	"""
	mode can be 'single', 'mosaic'
	msname  output measurementset
	nchan number of channel in spw
	centerfreq  centre frequency of spw
	chanwidth   channel width
	inttime      integration time 
	obstime   total time of observation
	numpointsource  Number of point sources to put in the field
	maxpointing error   each antenna will have get pointing error less than this value at each ppointing centre
	
	"""
	if(not os.path.exists(msname)):
		if(not createMS(msname=msname, nchan=nchan, centrefreq=centrefreq, chanwidth=chanwidth, inttime=inttime,obstime=obstime,noise=noise)):
			print(f"Could not create the measurement set {msname}")
			return
	u=simutil()
	# We make a copy to do baseline based source 
	tmpms=msname+"_temp.ms"
	tb=casatools.table()
	qa=casatools.quanta()
	sm=casatools.simulator()
	msmd=casatools.msmetadata()
	cl=casatools.componentlist()
	freq=qa.convert(qa.quantity(centrefreq), 'Hz')['value']
	dishdiam=17.0
	fieldofview=qa.quantity(consts.c.to('m/s').value/freq/dishdiam, 'rad')
	print(f"Fieldofview {fieldofview}")
	cmplist=f"Sources_{numpointsource}_{maxflux}.cl"
	msmd.open(msname)
	phcen=u.dir_m2s(msmd.phasecenter(0))
	msmd.done()
	
	minflux=f"{qa.div(qa.quantity(maxflux), 10)['value']}{qa.quantity(maxflux)['unit']}"
	createComponentList(clname=cmplist, refdir=phcen, maxflux=maxflux, minflux=minflux, numsource=numpointsource,
					 fov=qa.tos(fieldofview))
	bmEst=beamEstimator(msname=msname, field="0", spw="0")
	cl.open(cmplist)
	offsetdict=generateOffsetPerAntennaPerTime(msname=msname, maxoffset=maxpointError, distribution='uniform')
	for source in range(numpointsource):
		scl=casatools.componentlist()
		comp=cl.getcomponent(source)
		compdir=cl.getrefdir(source)
		comp={'component0':comp, 'nelements':1}
		scl.fromrecord(comp)
		shutil.rmtree("one_comp.cl", ignore_errors=True)
		scl.rename("one_comp.cl")
		scl.done()
		##now do this for everycomponent
		shutil.rmtree(tmpms, ignore_errors=True)
		shutil.copytree(msname, tmpms)
		tb.open(tmpms, nomodify=False)
		tb1=tb.taql("UPDATE "+tmpms+" SET DATA=0.0")
		tb1.done()
		tb1=tb.taql("UPDATE "+tmpms+" SET CORRECTED_DATA=0.0")
		tb1.done()
		tb.done()
		sm.openfromms(tmpms)
		sm.setdata(fieldid=[0], spwid=[0])
		sm.predict(complist='one_comp.cl')
		sm.done()
		dist=distFromCent(tmpms, field=0, compdirection=compdir)
		print(f"comp {source} dist={dist} ")
		tb.open(tmpms, nomodify=False)
		for ant_tim, offt in offsetdict.items():
			#time0=time.time()
			ant=ant_tim[0]
			tim=ant_tim[1]
			#distAntT=np.array(offt)+np.array(dist)
			distStr=[f"{offt[0]+dist[0]}rad", f"{offt[1]+dist[1]}rad"]
			voltFac=bmEst.estimateVoltageAtOffset(offset=distStr)
			tb1=tb.query(f"(ANTENNA1=={ant} || ANTENNA2 =={ant}) and TIME={tim}")
			dat=tb1.getcol('DATA')
			for chan in range(nchan):
				dat[:,chan,:] *= voltFac[chan]
			tb1.putcol('DATA',dat)
			#tb1.putcol("CORRECTED_DATA",dat)
			tb1.done()
			#print(f"1time for ant={ant}, t={tim}, {time.time()-time0}")
		tb1=tb.taql(f"UPDATE {msname}, {tmpms} t2 SET DATA=DATA+t2.DATA")
		tb1.done()
		tb1=tb.taql(f"UPDATE {msname}, {tmpms} t2 SET CORRECTED_DATA=CORRECTED_DATA+t2.DATA")
		tb1.done()
		tb.done()

def func(x,m,c):
  lsq=(m*x) + c
  return lsq 

def catalogue(image1,image2):
  file1=image1;file2=image2
  Img=bdsf.process_image(file1,  rms_box=(20,10))
  Img.write_catalog(catalog_type ='srl',format='fits')
  Img=bdsf.process_image(file2,  rms_box=(20,10))
  Img.write_catalog(catalog_type ='srl',format='fits')
  hdul1 = fits.open(file1[:-5]+'.pybdsf.srl.fits')
  data1 = hdul1[1].data
  Flux1= data1['Total_flux']
  hdul2 = fits.open(file2[:-5]+'.pybdsf.srl.fits')
  data2 = hdul2[1].data
  Flux2= data2['Total_flux']
  coeff, var = curve_fit(func,log10(Flux1),log10(Flux2))
  Flux1=array((Flux1))
  Flux2=array((Flux2))
  if Flux1.shape != Flux2.shape:
    # Determine which array is smaller and resize accordingly
    size_diff = len(Flux1) - len(Flux2)
    if size_diff > 0:
        # Truncate array1 to match array2's size
        Flux1 = Flux1[:len(Flux2)]
    else:
        # Truncate array2 to match array1's size
        Flux2 = Flux2[:len(Flux1)]
  flux_mean=mean((abs(Flux1-Flux2)/Flux1) * 100.0)
  variance = np.diagonal(var)
  SE = np.sqrt(variance)
  index=coeff[0]
  indexErr=SE[0]
  #spx.append(index)
  #spx_err.append(indexErr)
  yfit = func(log10(Flux1),coeff[0],coeff[1])
  fig = pl.figure(figsize=(10, 10))
  ax = pl.subplot(111)
  ax.plot(log10(Flux1),log10(Flux2),'o',ms=7,label=r'Differences % = ' + '%0.2f' % (flux_mean))
  ax.plot(log10(Flux1),yfit,'-',linewidth=2.0, label=r'$\alpha$'+' = ' '%0.2f $\pm$ %0.2f' % (index, indexErr))
  #matplotlib.pylab.autoscale(enable=True,tight=None)
  #pl.xlim(log10(Frequency[0]),log10(Frequency[2]))
  #ax.relim()
  ax.autoscale_view()
  box = ax.get_position()
  #ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
  #ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  ax.legend(numpoints=1,frameon=False)
  ax.tick_params(which='both',direction="in")
  ax.xaxis.set_minor_locator(AutoMinorLocator())
  ax.yaxis.set_minor_locator(AutoMinorLocator())
  #ax.autoscale()
  #fig.tight_layout()
  pl.xlabel('Image 1 flux (mJy)',size=14)
  pl.ylabel('Image 2 flux (mJy)',size=14)
  #pl.savefig('E_Relic_spectral_plot_1D.png')#,bbox_inches='tight')
  #pl.savefig('W_Relic_spectral_plot_1D.png')#,bbox_inches='tight')
  #pl.savefig('Galaxy_spectral_plot_1D.png')#,bbox_inches='tight')
  #pl.savefig('E_relic_flux_co-relation.png')#,bbox_inches='tight')
  #pl.savefig('W_relic_flux_co-relation.png')#,bbox_inches='tight')
  #pl.savefig('galaxy_flux_co-relation.png')#,bbox_inches='tight')
  pl.savefig(fits_dir+'/'+file1[:-5]+'_'+file2[:-5]+'_flux_comparisons.png')




def image_stat(image):
 inp_img=image
 fits_img=fits.open(inp_img)[0]
 image_data=fits_img.data
 image_header=fits_img.header
 max=np.max(image_data)
 sigclip = SigmaClip(sigma = 3.0)#, iters=5)
 image_data_sigclip = sigclip(image_data)
 image_min_clip, image_max_clip = np.min(image_data_sigclip),np.max(image_data_sigclip)
 image_mean_clip=np.mean(image_data_sigclip);image_median_clip=np.median(image_data_sigclip);image_std_clip=np.std(image_data_sigclip)
 image_mean, image_median, image_stddev = sigma_clipped_stats(image_data, sigma = 3.0)#, iters = 5)
 DR=max/image_stddev
 return DR,max,image_stddev

def syn_beam(image):
 inp_img=image
 fits_img=fits.open(inp_img)[0]
 image_data=fits_img.data
 image_header=fits_img.header
 return ((image_header['BMAJ']*3600), (image_header['BMIN']*3600), (image_header['BPA'] ))


def phase_corrupt1(msname: str="",deg: float="60.0"):
   tb=casatools.table()
   tb.open(msname)
   rad=float(deg)* (np.pi/180.)
   rad1=np.cos(rad)
   rad2=np.sin(rad)
   ###when antenna1 is 2 the multiply corrected_data with (cos60 +jsin60)

   tt=tb.taql('update '+str(msname)+' set DATA=iif(ANTENNA1==2,DATA*Complex('+str(rad1)+',' +str(rad2)+'), DATA) ')
   tt.done()

   ####when antenna2 is 2 multiply corrected_data with cos60-jsin60...note conjugate

   tt=tb.taql('update '+str(msname)+' set DATA=iif(ANTENNA2==2,DATA*Complex('+str(rad1)+',' +str(-rad2)+'), DATA) ')
   tt.done()
   tb.done()

def phase_corrupt2(msname: str="", deg: float="60.0"):
  tb=casatools.table()
  tb.open(msname, nomodify=False)
  rad=float(deg)* (np.pi/180.)
  rad1=np.cos(rad)
  rad2=np.sin(rad)
  ###get antenna1, antenna2 and corrected_data into python

  ant1=tb.getcol('ANTENNA1')
  ant2=tb.getcol('ANTENNA2')
  dat=tb.getcol('DATA')

  #### dat has shape (npol, nchan, nrow)

  ###for the rows where antenna1==2 multiply data with complex(cos60, sin60)

  dat[:,:, ant1==2] *= complex(rad1, rad2)

  ####for the rows where antenna2==2  multiply data with conjugate i.e complex(cos60, -sin60)
  dat[:,:, ant2==2] *= complex(rad1, -rad2)

  ###put modified data back

  tb.putcol('DATA', dat)

  tb.done()   

def phase_corrupt3(msname: str="", deg: float=""):
 tb.open(msname, nomodify=False)
###get antenna1, antenna2 and corrected_data into python
 ant1=tb.getcol('ANTENNA1')
 ant2=tb.getcol('ANTENNA2')
 dat=tb.getcol('DATA')
 nant=np.max(ant2)+1
#### dat has shape (npol, nchan, nrow)
 for i in range(nant):
#phase error in radian  radom between 20 and 60 degrees
    x=-deg;y=deg
    if x==-0 and y==0:
       phserr=0
    else:
       print(f"Phases between {x} and {y}")
       phserr=random.randrange(x,y)/180.0*np.pi
    #phserr=random.randrange(30,60)/180.0*np.pi
    phserr=deg*(np.pi/180.0)
    #phserr=0
    print(f"Applying phase {phserr} to antenna{i}")
    cphs=complex(np.cos(phserr), np.sin(phserr))
    ###for the rows where antenna1==1 multiply data with complex(cos(phserr), sin(phserr))
    dat[:,:, ant1==i] *= cphs
    ####for the rows where antenna2==i  multiply data with conjugate i.e
    dat[:,:, ant2==i] *= np.conj(cphs)
###put modified data back
 tb.putcol('DATA', dat)
 tb.done()   

def single_ant_phase_corrupt(msname: str, interval_minutes: int = 5, phase_range: tuple = (-10, 10), antenna_id: int = 2):
    """
    Adds random phase errors to data at regular intervals for a specific antenna in a measurement set.

    Parameters:
        msname (str): Name of the measurement set.
        interval_minutes (int): Time interval in minutes for applying phase errors.
        phase_range (tuple): Range of random phase errors in degrees (min, max).
        antenna_id (int): Antenna ID to which the phase errors will be applied.
    """
    tb=casatools.table()
    tb.open(msname, nomodify=False)
    try:
        # Retrieve relevant columns
        ant1 = tb.getcol('ANTENNA1')
        ant2 = tb.getcol('ANTENNA2')
        dat = tb.getcol('DATA')
        times = tb.getcol('TIME')  # Time in seconds
        
        # Define the interval in seconds
        interval_seconds = interval_minutes * 60
        start_time = times[0]
        next_time = start_time + interval_seconds
        # Initialize random phase values
        random_deg = np.random.uniform(*phase_range)
        rad = random_deg * (np.pi / 180.0)  # Convert to radians
        rad1 = np.cos(rad)
        rad2 = np.sin(rad)
        # Iterate over rows in time intervals
        for i, time in enumerate(times):
            if time >= next_time:  # Move to the next interval
                start_time = next_time
                next_time = start_time + interval_seconds
                # Generate a new random phase error for the new interval
                random_deg = np.random.uniform(*phase_range)
                rad = random_deg * (np.pi / 180.0)  # Convert to radians
                rad1 = np.cos(rad)
                rad2 = np.sin(rad)
            # Apply phase corruption for rows within the current interval
            if start_time <= time < next_time:
                if ant1[i] == antenna_id:
                    dat[:, :, i] *= complex(rad1, rad2)
                if ant2[i] == antenna_id:
                    dat[:, :, i] *= complex(rad1, -rad2)
        # Put modified data back into the measurement set
        tb.putcol('DATA', dat)

    finally:
        tb.done()

def all_ant_phase_corrupt(msname: str = "", deg: tuple = (-10, 10), interval_minutes: int = 5, refant: int = 3):
    """
    Adds random phase errors to all antennas except antenna 3 in a measurement set,
    applying new random errors every specified time interval.

    Parameters:
        msname (str): Name of the measurement set.
        deg (float): Maximum phase error range in degrees.
        interval_minutes (int): Time interval in minutes for applying phase errors.
        refant (int): Reference antenna number (in gaincal). 
    """
    if not msname:
        raise ValueError("Measurement set name must be provided.")

    tb=casatools.table()
    tb.open(msname, nomodify=False)

    try:
        # Retrieve columns
        ant1 = tb.getcol('ANTENNA1')
        ant2 = tb.getcol('ANTENNA2')
        dat = tb.getcol('DATA')
        times = tb.getcol('TIME')  # Time in seconds since the epoch

        # Calculate the interval in seconds
        interval_seconds = interval_minutes * 60

        # Determine time intervals
        start_time = times[0]
        end_time = times[-1]
        intervals = np.arange(start_time, end_time, interval_seconds)

        # Process data for each time interval
        for start, end in zip(intervals[:-1], intervals[1:]):
            # Generate a random phase error for each antenna in this interval
            phase_errors = {
                i: complex(np.cos(np.deg2rad(random.uniform(*deg))),
                           np.sin(np.deg2rad(random.uniform(*deg))))
                for i in range(np.max(ant2) + 1) if i != refant
            }

            print(f"Applying random phase errors from {start:.2f} to {end:.2f} seconds")

            # Identify rows within the current interval
            mask = (times >= start) & (times < end)

            # Apply phase errors
            for antenna, phase_error in phase_errors.items():
                phserr2=phase_error*(180/np.pi)
                print(f"Applying phase {phserr2} deg to antenna {antenna}")
                # Apply to ANTENNA1 rows
                dat[:, :, mask & (ant1 == antenna)] *= phase_error
                # Apply to ANTENNA2 rows with conjugate phase error
                dat[:, :, mask & (ant2 == antenna)] *= np.conj(phase_error)

        # Write the modified data back
        tb.putcol('DATA', dat)

    finally:
        tb.done()


def all_ant_phase_amp_corrupt(msname: str = "", deg: tuple = (-10, 10), amp_frac: tuple = (0.01, 0.10),
                          interval_minutes: int = 5, refant: int = 3):
    """
    Adds random complex gain errors (phase and amplitude) to all antennas except a reference antenna,
    applying new random errors every specified time interval.

    Parameters:
        msname (str): Name of the measurement set.
        deg (tuple): Phase error range in degrees (min, max).
        amp_frac (tuple): Amplitude error range as a fraction (e.g., (0.01, 0.10) for 1%-10%).
        interval_minutes (int): Time interval in minutes for applying gain errors.
        refant (int): Reference antenna number to skip corruption.
    """
    if not msname:
        raise ValueError("Measurement set name must be provided.")

    tb = casatools.table()
    tb.open(msname, nomodify=False)

    try:
        # Get antenna and data info
        ant1 = tb.getcol('ANTENNA1')
        ant2 = tb.getcol('ANTENNA2')
        dat = tb.getcol('DATA')  # Shape: (nchan, ncor, nrow)
        times = tb.getcol('TIME')

        interval_seconds = interval_minutes * 60
        start_time, end_time = times[0], times[-1]
        intervals = np.arange(start_time, end_time, interval_seconds)

        nant = np.max([ant1.max(), ant2.max()]) + 1

        for start, end in zip(intervals[:-1], intervals[1:]):
            # Generate random complex gains per antenna
            gain_errors = {}
            for ant in range(nant):
                if ant == refant:
                    continue
                phase_deg = random.uniform(*deg)
                amp_error = 1.0 + random.uniform(-amp_frac[1], amp_frac[1])  # e.g., 1.07 or 0.95
                gain = amp_error * np.exp(1j * np.deg2rad(phase_deg))
                gain_errors[ant] = gain

            print(f"\nApplying errors from {start:.2f} to {end:.2f} seconds")

            # Apply to all data rows within time interval
            mask = (times >= start) & (times < end)

            for ant, gain in gain_errors.items():
                amp_percent = (abs(gain) - 1.0) * 100
                phase_deg = np.angle(gain, deg=True)
                print(f"  Antenna {ant}: phase = {phase_deg:+.2f} deg, amp = {amp_percent:+.1f}%")

                dat[:, :, mask & (ant1 == ant)] *= gain
                dat[:, :, mask & (ant2 == ant)] *= np.conj(gain)

        # Write modified data back
        tb.putcol('DATA', dat)

    finally:
        tb.done()

def channel_phase_corrupt(msname: str = "", deg: tuple = (-10, 10), interval_minutes: int = 5):
    """
    Adds random phase errors to each frequency channel in a measurement set,
    applying new random errors every specified time interval.

    Parameters:
        msname (str): Name of the measurement set.
        deg (tuple): Phase error range in degrees (min, max).
        interval_minutes (int): Time interval in minutes for applying new phase errors.
    """
    if not msname:
        raise ValueError("Measurement set name must be provided.")

    tb = casatools.table()
    tb.open(msname, nomodify=False)

    try:
        # Retrieve columns
        dat = tb.getcol('DATA')   # Shape: (nchan, ncor, nrow)
        times = tb.getcol('TIME')  # Shape: (nrow,)

        nchan = dat.shape[0]

        # Convert interval to seconds
        interval_seconds = interval_minutes * 60

        # Define time bins
        start_time = times.min()
        end_time = times.max()
        time_bins = np.arange(start_time, end_time, interval_seconds)

        for t_start, t_end in zip(time_bins[:-1], time_bins[1:]):
            print(f"Applying phase error to channels between {t_start:.2f} and {t_end:.2f} seconds")

            # Get rows in this time interval
            mask = (times >= t_start) & (times < t_end)

            # Generate random phase errors for each channel
            phase_errors = np.array([
                np.exp(1j * np.deg2rad(random.uniform(*deg)))
                for _ in range(nchan)
            ])

            for ch in range(nchan):
                # Apply the phase error to all polarizations and all rows in this interval
                dat[ch, :, mask] *= phase_errors[ch]

                # For debugging
                ph_deg = np.angle(phase_errors[ch], deg=True)
                print(f"  Channel {ch}: phase {ph_deg:.2f}°")

        # Write modified data back
        tb.putcol('DATA', dat)

    finally:
        tb.done()

def channel_amp_phase_corrupt(
    msname: str = "",
    deg: tuple = (-10, 10),
    amp_frac: tuple = (0.01, 0.10),
    interval_minutes: int = 5
 ):
    """
    Adds random amplitude and phase errors per frequency channel to a Measurement Set,
    applying new random errors for each time interval.

    Parameters:
        msname (str): Name of the measurement set.
        deg (tuple): Phase error range in degrees (min, max).
        amp_frac (tuple): Amplitude error fraction range (e.g., 0.01 to 0.10 for 1%-10%).
        interval_minutes (int): Time interval in minutes to update channel-wise gain errors.
    """
    if not msname:
        raise ValueError("Measurement set name must be provided.")

    tb = casatools.table()
    tb.open(msname, nomodify=False)

    try:
        dat = tb.getcol('DATA')       # shape: (nchan, npol, nrow)
        times = tb.getcol('TIME')     # shape: (nrow,)
        nchan = dat.shape[0]

        # Calculate time intervals
        interval_seconds = interval_minutes * 60
        start_time, end_time = times.min(), times.max()
        time_bins = np.arange(start_time, end_time, interval_seconds)

        for t_start, t_end in zip(time_bins[:-1], time_bins[1:]):
            mask = (times >= t_start) & (times < t_end)
            if not np.any(mask):
                continue

            print(f"\nApplying per-channel corruption between {t_start:.2f} and {t_end:.2f} seconds")

            # Create per-channel gain factors (complex)
            gains = []
            for ch in range(nchan):
                phase_deg = random.uniform(*deg)
                amp_scale = 1.0 + random.uniform(-amp_frac[1], amp_frac[1])
                gain = amp_scale * np.exp(1j * np.deg2rad(phase_deg))
                gains.append(gain)
                print(f"  Channel {ch:3d}: phase = {phase_deg:+6.2f}°, amp = {((amp_scale-1)*100):+5.1f}%")

            gains = np.array(gains, dtype=complex).reshape((nchan, 1, 1))  # shape: (nchan, 1, 1)
            dat[:, :, mask] *= gains  # Broadcasting over pol and row

        # Write modified data back
        tb.putcol('DATA', dat)

    finally:
        tb.done()

def phase_amp_corrupt(msname: str="", deg: float="60.0"):
 tb.open(msname, nomodify=False)
###get antenna1, antenna2 and corrected_data into python
 ant1=tb.getcol('ANTENNA1')
 ant2=tb.getcol('ANTENNA2')
 dat=tb.getcol('DATA')
 nant=np.max(ant2)+1
#### dat has shape (npol, nchan, nrow)
 for i in range(nant):
#phase error in radian  radom between 20 and 60 degrees
    x=deg-5;y=deg+5
    phserr=random.randrange(x,y)/180.0*np.pi
    #phserr=random.randrange(30,60)/180.0*np.pi
    #phserr=deg*(np.pi/180.0)
    #phserr=0
    print(f"Applying phase {phserr} to antenna{i}")
    cphs=complex(np.cos(phserr), np.sin(phserr))*random.uniform(0.1,0.10)
    ###for the rows where antenna1==1 multiply data with complex(cos(phserr), sin(phserr))
    dat[:,:, ant1==i] *= cphs
    ####for the rows where antenna2==i  multiply data with conjugate i.e
    dat[:,:, ant2==i] *= np.conj(cphs)
###put modified data back
 tb.putcol('DATA', dat)
 tb.done()   
		

def antenna_list(vis):

  ms_file = vis
  #tb = table()
  antenna_table = ms_file + "/ANTENNA"
  tb.open(antenna_table)
  antenna_names = tb.getcol("NAME")
  antenna_positions = tb.getcol("POSITION")
  tb.close()
  center_position = np.mean(antenna_positions, axis=1)
  distances = np.linalg.norm(antenna_positions - center_position[:, np.newaxis], axis=0)
  sorted_indices = np.argsort(distances)
  sorted_antenna_names = antenna_names[sorted_indices]
  sorted_antenna_names = ','.join(sorted_antenna_names) 
  
  return sorted_antenna_names

def antenna_list2(vis):
    """
    Returns a comma-separated list of antenna numbers sorted by their distance
    from the center of the array, along with their sorted positions.

    Parameters:
    vis (str): Path to the Measurement Set (MS) file.

    Returns:
    str: Comma-separated sorted antenna numbers.
    """
    try:
        # Initialize CASA table tool
        #tb = table()

        # Open the ANTENNA subtable
        antenna_table = f"{vis}/ANTENNA"
        tb.open(antenna_table)

        # Get antenna names and positions
        antenna_names = tb.getcol("NAME")
        antenna_positions = tb.getcol("POSITION")

        # Close the table
        tb.close()

        # Calculate the center position
        center_position = np.mean(antenna_positions, axis=1)

        # Compute distances of antennas from the center position
        distances = np.linalg.norm(antenna_positions - center_position[:, np.newaxis], axis=0)

        # Sort antenna indices by their distances
        sorted_indices = np.argsort(distances)

        # Generate antenna numbers (starting from 1)
        antenna_numbers = np.arange(1, len(antenna_names) + 1)

        # Sort the antenna numbers based on distances
        sorted_antenna_numbers = antenna_numbers[sorted_indices]

        # Convert the sorted list of numbers to a comma-separated string
        return ','.join(map(str, sorted_antenna_numbers))
    
    except Exception as e:
        print(f"Error processing antenna list: {e}")
        return ""
  


def baseline(ant_file,frequency):
 data=np.loadtxt(ant_file,usecols=(0,1))
 x=data[:,0];y=data[:,1]
 antenna_positions = list(zip(x, y))
 
 def euclidean_distance(p1, p2):
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
        
 max_distance = 0
 for i in range(len(antenna_positions)):
   for j in range(i + 1, len(antenna_positions)):
       dist = euclidean_distance(antenna_positions[i], antenna_positions[j])
       if dist > max_distance:
            max_distance = dist 
 long_dist=max_distance
 #print(long_dist)
 c=3e8 #in cm ;3e8 in m
 freq=float(frequency)
 diameter=long_dist
 wave=c/freq
 resolution = wave/diameter * (180/np.pi)
 resolution_asec=resolution*3600
 return resolution_asec

def VLA_sensitivity(npol=2,N=25,tint=3600,v=1e9):
 nc=0.92
 npol=npol
 N=N
 tint=tint*3600
 v=v
 SEFD=420
 S=SEFD/ (nc * sqrt(npol * N*(N-1) *tint *v) )
 return S

def simulation_noise(vis,t_track):
   msmd.open(vis)
   N_ant=msmd.nantennas()
   nbaseline=N_ant*(N_ant-1)/2.
   tint=msmd.effexposuretime()
   tint=tint['value']
   msmd.done()
   v=msmd.bandwidths()[0]
   npol=msmd.ncorrforpol()[0]
   N=nchan = msmd.nchan()
   ms.open(vis)
   spwInfo = ms.getspectralwindowinfo()
   print(spwInfo)
   freq=spwInfo['0']['Chan1Freq']
   chan_width=spwInfo['0']['ChanWidth']
   nchan=spwInfo['0']['NumChan']
   sensitivity=VLA_sensitivity(npol=npol,N=nchan,tint=tint,v=v)						 
   msmd.done()
   t_track=t_track #3600sec
   scale_factor=0.83/sqrt(t_track/1.0)
   rms_exp=scale_factor*sqrt(npol*nbaseline*nchn*nint) #nint=total on-source time / integration time=3600/5
   return rms_exp

def cleaning(vis,imagename,datacolumn,imagesize,cellsize,fits_dir):
 clean_output=tclean(vis=vis, imagename=imagename, datacolumn=datacolumn,weighting='briggs',robust=0,
      imsize=imagesize, cell=cellsize, pblimit=-0.01, niter=10000, savemodel='modelcolumn')   
 fits_output=exportfits(imagename=imagename+'.image',fitsimage=fits_dir+'/'+str(imagename)+'.image'+'.fits',overwrite=True) 
 return clean_output,fits_output
 
  
def cleaning1(vis,imagename,column):
  clean_output=tclean(vis=vis, imagename=imagename,field=img_fields,spw=spw,
       gridder='standard',wprojplanes=1, pblimit=-0.1, imsize=imagesize, cell=cellsize, specmode='mfs',
       deconvolver='hogbom', nterms=nterms, scales=scale, datacolumn=column, smallscalebias=0.9, phasecenter=phasecenter, restoration=True, restoringbeam='common', 
       interactive=False, niter=1000,  weighting='briggs',robust=0.0,
       stokes='I',threshold='1e-6Jy',sidelobethreshold=2.0,nsigma=4.0,cyclefactor=3.0,parallel=parallel,calcpsf=True, calcres=True,restart=True)
  fits_out=exportfits(imagename=imagename+'.image',fitsimage=fits_dir+'/'+imagename+'.fits',overwrite=True)

  return clean_output,fits_out
  
def cleaning2(vis,imagename,column):
  clean_output=tclean(vis=vis, imagename=imagename+'_wproj',field=img_fields,spw=spw,
       gridder='wproject',wprojplanes=1, pblimit=-0.1, imsize=imagesize, cell=cellsize, specmode='mfs',
       deconvolver='hogbom', nterms=nterms, scales=scale, datacolumn=column, smallscalebias=0.9, phasecenter=phasecenter, restoration=True, restoringbeam='common', 
       interactive=False, niter=1000,  weighting='briggs',robust=0.0,
       stokes='I', threshold='1e-6Jy',sidelobethreshold=2.0,nsigma=4.0,cyclefactor=3.0,parallel=parallel,calcpsf=True, calcres=True,restart=True)
  fits_out=exportfits(imagename=imagename+'_wproj'+'.image',fitsimage=fits_dir+'/'+imagename+'_wproj'+'.fits',overwrite=True)

  return clean_output,fits_out


def cleaning3(vis,imagename,column):
  clean_output=tclean(vis=vis, imagename=imagename+'_awp2',field=img_fields,spw=spw,
       gridder='awp2',wprojplanes=1, pblimit=-0.1, imsize=imagesize, cell=cellsize, specmode='mfs',
       deconvolver='hogbom', nterms=nterms, scales=scale, datacolumn=column, smallscalebias=0.9, phasecenter=phasecenter, restoration=True, restoringbeam='common', 
       interactive=False, niter=1000,  weighting='briggs',robust=0.0,
       stokes='I', threshold='1e-6Jy',sidelobethreshold=2.0,nsigma=4.0,cyclefactor=3.0,parallel=parallel,calcpsf=True, calcres=True,restart=True)
  fits_out=exportfits(imagename=imagename+'_awp2'+'.image',fitsimage=fits_dir+'/'+imagename+'_awp2'+'.fits',overwrite=True)

  return clean_output,fits_out

def selfcal(vis,caltable,solint): 
 gaincal(vis=vis,
        caltable=caltable, spw='',
        solint=solint, combine='',field='',selectdata=True,solnorm=False,
        refant='3',   
        minblperant=4, minsnr=5.0, gaintype='G', calmode='p', append=False,
        parang=False,
        gaintable=[])  
 plotms(vis=caltable,xaxis='Time',yaxis='phase',iteraxis='antenna',gridrows=2,gridcols=2,
         plotfile=caltable+'_phase_error.png')       
 applycal(vis=vis,spw='',selectdata=True,
         gaintable=[caltable], field='',interp=['linear'], calwt=False, parang=False,
         applymode='calonly')

def simulation(antfile,project_name,model,mapsize):
   sim_output=simobserve(project = project_name,
           skymodel = model,
           inbright='',
           #indirection='J2000 359.9994708deg 30.0009052deg',  #J2000 19h00m00 -40d00m00
           incell='',
           incenter='',
           inwidth='100MHz',
           setpointings = True,
           #ptgfile='',
           #direction="J2000 359.9994708deg 30.0009052deg",
           refdate="2023/02/03",
           mapsize = mapsize,
           integration = '60s',
           obsmode = 'int',
           antennalist = antfile,
           hourangle = 'transit',
           totaltime = '1800s',
           thermalnoise = '',
           outframe='LSRK',
           graphics = 'none')  
           
   return sim_output
   

def plotting(vis):
 plot_out1=plotms(vis=vis,xaxis="uvdist",yaxis="Amp",showgui=False,ydatacolumn="data",plotfile='plot_files/'+vis+'_uvdist_vs_ampplot.png',clearplots=True,overwrite=True,avgchannel='0',avgtime='0')
 plot_out2=plotms(vis=vis,xaxis="Time",yaxis="Amp",showgui=False,ydatacolumn="data",plotfile='plot_files/'+vis+'_time_amp_corrmodel.png',clearplots=True,overwrite=True,avgchannel='0',avgtime='0')  
 plot_out3=plotms(vis=vis,xaxis="Channel",yaxis="Amp",showgui=False,ydatacolumn="data",plotfile='plot_files/'+vis+'_chn_amp_corrmodel.png',clearplots=True,overwrite=True,avgchannel='0',avgtime='0') 
 plot_out4=plotms(vis=vis,xaxis="U",yaxis="V",showgui=False,ydatacolumn="data",plotfile='plot_files/'+vis+'_u_v_corrmodel.png',clearplots=True,overwrite=True,avgchannel='0',avgtime='0') 
 return plot_out1,plot_out2,plot_out3,plot_out4
 
 
def manual_simulation(conf_file,configdir,vis,freq,deltafreq,freqresolution,nchannels,integrationtime,total_obs_time,model_img_name,RA,DEC):
  conf_file = conf_file
  configdir = configdir
  u = simutil()
  xx,yy,zz,diam,padnames,padnames2,telescope,posobs = u.readantenna(configdir+conf_file)
  ms_name = vis    ## Name of your measurement set
  sm.open(ms_name)
  pos_ngVLA = me.observatory('ngvla')
  sm.setconfig(telescopename = telescope, x = xx, y = yy, z = zz,
                    dishdiameter = diam.tolist(), mount = 'alt-az',
                    antname = padnames, padname = padnames,
                    coordsystem = 'global', referencelocation = pos_ngVLA)
  sm.setspwindow(spwname = 'Band1', freq = freq, deltafreq = deltafreq,
  freqresolution = freqresolution, nchannels = nchannels, stokes = 'RR RL LR LL')
  sm.setfeed('perfect R L')
  sm.setfield(sourcename = 'My source',
                  sourcedirection = ['J2000',str(RA)+'deg',str(DEC)+'deg'])
  sm.setlimits(shadowlimit = 0.001, elevationlimit = '8.0deg')
  sm.setauto(autocorrwt = 0.0)                
  integrationtime = integrationtime
  sm.settimes(integrationtime = integrationtime, usehourangle = True,
            referencetime = me.epoch('utc', 'today'))
  starttime = '0h'
  stop_time = float(total_obs_time)/3600.0
  stoptime=str(stop_time)+'h'
  sm.observe('My source', 'Band1', starttime = starttime, stoptime = stoptime)
  importfits(fitsimage = model_img_name, imagename = model_img_name+'.image',overwrite=True)    
  sm.openfromms(vis)
  sm.predict(imagename = model_img_name+'.image') 
  #sigma_simple = '1.4mJy'
  #os.system('cp -r ngVLA_214_ant_60s_noise_free.ms ngVLA_214_ant_60s_noisy.ms')
  sm.openfromms(ms_name)
  #sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)
  #sm.corrupt()
  sm.done()
  sm.close()
  #return vis
