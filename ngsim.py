import sys,os
sys.path.append(os.getcwd())

from ngsim_helpers import *
configdir ='/home/rarg/software/casa-6.6.4-34-py3.8.el8/data/alma/simmos/'

Simulation=False
Plotting=False
Vis_calculations=False
noise_corruption=False
phase_corruption=False
phaseamp_corrupt=False
chan_phase_corruption=False
chan_amp_phase_corruption=False
Cleaning=True
Cleaning1=False
Cleaning2=False
Cleaning3=False
directories_rm=False

#simulation parameters
project_name='ngVLA_HDR'
centre_freq='2e9'    #in Hz
total_channels=5 
total_bw='1e9'       #in Hz
total_obs_time='3600'   #in sec
deltafreq = float(total_bw)/total_channels
freq=float(centre_freq) - (float(total_bw) / 2.) + (float(deltafreq) / 2.)
freq=str(freq)
freqresolution = str(deltafreq)
deltafreq= str(deltafreq)
integrationtime='10s'
inp_img='EGSPS8.4.FITS'
observatory='ngVLA'
conf_file='ngvla-revD.main.cfg'
antfile=['ngvla-revD.main.cfg']
#antfile=(('ngvla-revD.core.cfg','ngvla-revD.mid.cfg','ngvla-revD.spiral_core.cfg','ngvla-revD.lba.cfg','ngvla-revD.sba.cfg','ngvla-revD.main.cfg','ngvla-revD.spiral.cfg'))
string='myrun_'
my_string='model_corr_'
point_err=0 
interval_minutes=5
phase_err=[5]
number=i=0 #0=5min,1=10min, 2=15min, 3=20min, 4=30min

#model parameters
fits_img=fits.open(inp_img)[0]
data=fits_img.data
Header=fits_img.header
if len(data.shape)==4:
 (a,b,XX,YY)=data.shape
if len(data.shape)==3:
 (a,XX,YY)=data.shape
if len(data.shape)==2:
 (XX,YY)=data.shape
delta1=abs(Header['CDELT1'])
RA=Header['CRVAL1']
#DEC=[60,40,20,10,0,-10,-20,-30,-40]
DEC=Header['CRVAL2']


#imaging parameters
imagesize=720
cellsize='0.22arcsec'
img_fields='0'
spw=''       
nterms=2
scale=[0,5,12]
#phasecenter='J2000 17:45:28.0297 30.00.00'
#phasecenter='J2000 '+str(ra)+'deg '+str(dec)+'deg'
phasecenter=''
parallel=False

for k in ((phase_err)):
 for i in ((antfile)):
   phase_err=k
   model=inp_img[:-5]+'_'+str(RA)+'deg'+'_'+str(DEC)+'deg'+'.fits'
   print(model)
   ant_file=i
   delta=baseline(ant_file,centre_freq)       
   header(inp_img,RA,DEC,delta)
   model_img_name=inp_img[:-5]+'_'+str(RA)+'deg'+'_'+str(DEC)+'deg'+'.fits'
   project_name=inp_img[:-5]+'_'+ant_file[:-4]+'_'+str(point_err)+'_'+str(phase_err)+'_'+str(number)
   imagename=inp_img+'_'+ant_file[:-4]+'_clean_image'
   if not os.path.exists(project_name):
      os.mkdir(str(project_name))
   vis=str(project_name)+'/'+str(project_name)+'.'+ant_file[:-4]+'.ms'
   fits_dir=model[:-5]+'_simulation_'+str(observatory)+'_'+str(point_err)+'_'+str(phase_err)+'deg_'+str(my_string)+str(number)
   print("visibility name = "+str(vis))
   
   if Simulation==True: 
     #simulation(ant_file,project_name,model,mapsize)
     manual_simulation(conf_file,configdir,vis,freq,deltafreq,freqresolution,total_channels,integrationtime,total_obs_time,model_img_name,RA,DEC) 
  
   if noise_corruption==True:
     sm.openfromms(vis)
     sigma_simple = '2.5mJy'
     sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)
     sm.corrupt()
     sm.done()						 
  
   if phase_corruption==True:
     #phase_corrupt3(msname=vis,deg=phase_err)
     #print("phase_corruption is True")
     #print(vis)	
     all_ant_phase_corrupt(msname=vis, deg=(-phase_err,phase_err), interval_minutes=interval_minutes,refant=3)       
  
   if phaseamp_corrupt==True:
     #phase_amp_corrupt(msname=vis,deg=phase_err)
     all_ant_phase_amp_corrupt(msname=vis, deg=(-phase_err,phase_err), amp_frac=(0.01, 0.10), interval_minutes=interval_minutes,refant=3) 
   
   if chan_phase_corruption==True:
     channel_phase_corrupt(msname=vis, deg= (-phase_err,phase_err), interval_minutes=interval_minutes)
   
   if chan_amp_phase_corruption==True:
    channel_amp_phase_corrupt(msname=vis,deg= (-phase_err,phase_err),amp_frac=(0.01, 0.10),interval_minutes = interval_minutes)
     
   if Vis_calculations==True:
    vis=str(project_name)+'/'+str(project_name)+'.'+ant_file[:-4]+'.ms'
    psfsize=psf_size(vis)
    #print(psfsize[0])
    print("The resolution at given frequency is "+str(psfsize[1])+' asec')
    #print (psfsize[1]/3.)
    cellsize=psfsize[1]/4.
    print ("The required cellsize is "+str(cellsize)+' asec')
    image_size=fov(vis)
    print("The FOV at given frequency is "+str(image_size)+' amin')
    imsize=(4*image_size*60)/cellsize
    print("The required image size is "+str(imsize))
    N_wplane=((4*image_size*60)/psfsize[1]) * ((4*image_size*60)/206265.0)  #1radian=206265"
    print ("The required w-planes are "+str(N_wplane))
  
   if not os.path.exists(fits_dir):
     os.makedirs(fits_dir)
  
   if Plotting==True:
    plotting(vis=vis)
  
   if Cleaning==True:
     fd=open(fits_dir+"/"+"image_statistics.txt", "w")
     antennas=antenna_list2(vis)
     cellsize='0.22arcsec'
     cleaning(vis,imagename,'data',imagesize,cellsize,fits_dir)  
     solints=['inf','15min','5min','int']
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.selfcal'
      selfcal(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      imagename=imagename+'.sc.'+str(solint)
      cleaning(vis,imagename,'corrected',imagesize,cellsize,fits_dir) 
      stat=image_stat(fits_dir+'/'+str(imagename)+'.image'+'.fits')
      beam=syn_beam(fits_dir+'/'+str(imagename)+'.image'+'.fits')
      fd.write(f"Image: {imagename+'.image'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
      fd.close()
  
   if Cleaning1==True:
    cleaning1(vis=vis,imagename=imagename,column='data')
    selfcal(vis) 
    cleaning1(vis=vis,imagename=imagename+'_sc',column='corrected')
    shutil.rmtree(vis+'.G.selfcal')
  
   if Cleaning2==True:
    cleaning2(vis=vis,imagename=imagename,column='data')
    selfcal(vis) 
    cleaning2(vis=vis,imagename=imagename+'_sc',column='corrected')
    shutil.rmtree(vis+'.G.selfcal')
   
   directory = vis
   if os.path.isdir(directory):
     print("Yes, the directory exists.")
   else:
     print("No, the directory does not exist.")

   if directories_rm==True:
     #shutil.rmtree(vis);shutil.rmtree(project_name)
     directories=glob.glob(imagename+'**')
     for directory in directories:
      if os.path.isdir(directory):
           shutil.rmtree(directory)   
