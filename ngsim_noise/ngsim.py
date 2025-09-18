import sys,os
sys.path.append(os.getcwd())

from ngsim_helpers import *

Simulation=False
Plotting=False
Vis_calculations=False
noise_corruption=True
phase_corruption=True
phaseamp_corrupt=False
chan_phase_corruption=True
chan_amp_phase_corruption=False
Cleaning=False
Cleaning1=True
Cleaning1_sc_ap=False
Cleaning1_cube=False
Cleaning1_sc_ap_cube=False
Cleaning2=False
Cleaning3=False
directories_rm=True

#simulation parameters
project_name='ngVLA_DR'
centre_freq='2e9'    #in Hz
total_channels=10
total_bw='1e9'       #in Hz
total_obs_time='21600s'   #in sec
deltafreq = float(total_bw)/total_channels #channel_width
freq=float(centre_freq) - (float(total_bw) / 2.) + (float(deltafreq) / 2.)
freq_end = freq + (total_channels - 1) * deltafreq
freq=str(freq)
freqresolution = str(deltafreq)
deltafreq= str(deltafreq)
integrationtime='10s'
inp_img='EGSPS8.4.FITS'
observatory='ngVLA'
conf_file='ngvla-revD.main.cfg'
antfile=['ngvla-revD.main.cfg']
#antfile=(('ngvla-revD.core.cfg','ngvla-revD.mid.cfg','ngvla-revD.spiral_core.cfg','ngvla-revD.lba.cfg','ngvla-revD.sba.cfg','ngvla-revD.main.cfg','ngvla-revD.spiral.cfg'))
string='myrun_noise_'
my_string='noise_'
point_err=0 
interval_min_error=[5]#,10,15,30]
phase_error=[5]#,10,20,25,40]
number=i=1 #0=5min,1=10min, 2=15min, 3=20min, 4=30min



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
alpha=-0.7

#imaging parameters
imagesize=1000
cellsize='0.22arcsec'
img_fields='0'
spw=''       
nterms=2
scale=[0,1,2]
#phasecenter='J2000 17:45:28.0297 30.00.00'
#phasecenter='J2000 '+str(ra)+'deg '+str(dec)+'deg'
phasecenter=''
parallel=False
#solints=['inf','30min','5min','int']
solints_phase=['inf','30min','5min','int']
solints_ap=['int']

for t in interval_min_error:
 for k in phase_error:
  for i in antfile:
   phase_err=k;interval_minutes=t;number=t
   model=inp_img[:-5]+'_'+str(RA)+'deg'+'_'+str(DEC)+'deg_cube'+'.fits'
   print(model)
   ant_file=i
   delta=baseline(ant_file,centre_freq)       
   #delta=delta1*3600
   #header(inp_img,RA,DEC,delta)
   cube_generate(inp_img,RA,DEC,delta,alpha,centre_freq,freq,freq_end,total_channels)
   model_img_name=inp_img[:-5]+'_'+str(RA)+'deg'+'_'+str(DEC)+'deg_cube'+'.fits'
   project_name=inp_img[:-5]+'_'+ant_file[:-4]+'_'+str(point_err)+'_'+str(phase_err)+'_'+str(my_string)+str(number)
   imagename=inp_img+'_'+ant_file[:-4]+'_clean_image'
   if not os.path.exists(project_name):
      os.mkdir(str(project_name))
   vis=str(project_name)+'/'+str(project_name)+'.'+ant_file[:-4]+'.ms'
   fits_dir=model[:-5]+'_simulation_'+str(observatory)+'_'+str(point_err)+'_'+str(phase_err)+'deg_'+str(my_string)+str(number)
   fits_dir2=model[:-5]+'_simulation_'+str(observatory)+'_'+str(point_err)+'_'+str(phase_err)+'deg_'+str(my_string)+str(number)+'_img_cube'
   print("visibility name = "+str(vis))
   
   if Simulation==True: 
     mapsize='' #f'{delta1*3600*XX}arcsec'
     simulation(ant_file,project_name,model,mapsize,integrationtime,total_obs_time)
     #manual_simulation(conf_file,vis,freq,deltafreq,freqresolution,total_channels,integrationtime,total_obs_time,model_img_name,RA,DEC) 
  
   if noise_corruption==True:
     t_track=int(total_obs_time.rstrip('s'))/3600.0
     sensitivity=0.44
     sigma_simple=simulation_noise(vis,t_track,sensitivity)
     sigma_simple=str(sigma_simple)+'mJy'
     print(sigma_simple)
     sm.openfromms(vis)
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
   
   if not os.path.exists(fits_dir2):
     os.makedirs(fits_dir2)
   
   if Plotting==True:
    plotting(vis=vis)
  
   if Cleaning==True:
     fd=open(fits_dir+"/"+"image_statistics.txt", "w")
     antennas=antenna_list2(vis)
     cellsize='0.22arcsec'
     if os.path.exists(imagename+'.mask'):
        os.system("rm -rf "+str(imagename+'.mask')) 
     cleaning(vis,imagename,'data',imagesize,cellsize,fits_dir)  
     stat=image_stat(fits_dir+'/'+str(imagename)+'.image'+'.fits')
     beam=syn_beam(fits_dir+'/'+str(imagename)+'.image'+'.fits')
     fd.write(f"Image: {imagename+'.image'+'.fits'}\n")  # Write the image name
     fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
     fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")

     solints=solints_phase
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.selfcal'
      selfcal(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      sc_imagename=imagename+'.sc.'+str(solint)
      if os.path.exists(sc_imagename+'.mask'):
        os.system("rm -rf "+str(sc_imagename+'.mask')) 
      cleaning(vis,sc_imagename,'corrected',imagesize,cellsize,fits_dir) 
      stat=image_stat(fits_dir+'/'+str(sc_imagename)+'.image'+'.fits')
      beam=syn_beam(fits_dir+'/'+str(sc_imagename)+'.image'+'.fits')
      fd.write(f"Image: {sc_imagename+'.image'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
     fd.close()
  
   if Cleaning1==True:
     fd=open(fits_dir+"/"+"image_statistics.txt", "w")
     antennas=antenna_list2(vis)
     cellsize='0.22arcsec'
     if os.path.exists(imagename+'.mask'):
        os.system("rm -rf "+str(imagename+'.mask')) 
     cleaning1(vis,imagename,'data',imagesize,cellsize,fits_dir,scale,nterms,parallel)  
     stat=image_stat(fits_dir+'/'+str(imagename)+'.image.tt0'+'.fits')
     beam=syn_beam(fits_dir+'/'+str(imagename)+'.image.tt0'+'.fits')
     fd.write(f"Image: {imagename+'.image.tt0.'+'.fits'}\n")  # Write the image name
     fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
     fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")

     solints=solints_phase
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.selfcal'
      selfcal(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      sc_imagename=imagename+'.sc.'+str(solint)
      if os.path.exists(sc_imagename+'.mask'):
        os.system("rm -rf "+str(sc_imagename+'.mask')) 
      cleaning1(vis,sc_imagename,'corrected',imagesize,cellsize,fits_dir,scale,nterms,parallel)
      stat=image_stat(fits_dir+'/'+str(sc_imagename)+'.image.tt0'+'.fits')
      beam=syn_beam(fits_dir+'/'+str(sc_imagename)+'.image.tt0'+'.fits')
      fd.write(f"Image: {sc_imagename+'.image.tt0.'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
     fd.close()
  
   if Cleaning1_sc_ap==True:
     fd=open(fits_dir+"/"+"image_statistics.txt", "a")
     cellsize='0.22arcsec'
     solints=solints_ap
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.ap.selfcal'
      selfcal_ap(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      sc_imagename=imagename+'.sc.ap.'+str(solint)
      if os.path.exists(sc_imagename+'.mask'):
        os.system("rm -rf "+str(sc_imagename+'.mask')) 
      cleaning1(vis,sc_imagename,'corrected',imagesize,cellsize,fits_dir,scale,nterms,parallel)
      stat=image_stat(fits_dir+'/'+str(sc_imagename)+'.image.tt0'+'.fits')
      beam=syn_beam(fits_dir+'/'+str(sc_imagename)+'.image.tt0'+'.fits')
      fd.write(f"Image: {sc_imagename+'.image.tt0'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
     fd.close()

   if Cleaning1_cube==True:
     fd=open(fits_dir2+"/"+"image_statistics.txt", "w")
     antennas=antenna_list2(vis)
     cellsize='0.22arcsec'
     if os.path.exists(imagename+'.mask'):
        os.system("rm -rf "+str(imagename+'.mask')) 
     cleaning1_cube(vis,imagename,'data',imagesize,cellsize,fits_dir2,scale,nterms,parallel,total_channels,freq,deltafreq)  
     stat=image_stat(fits_dir2+'/'+str(imagename)+'.image'+'.fits')
     beam=syn_beam(fits_dir2+'/'+str(imagename)+'.image'+'.fits')
     fd.write(f"Image: {imagename+'.image'+'.fits'}\n")  # Write the image name
     fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
     fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")

     solints=solints_phase
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.selfcal'
      selfcal(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      sc_imagename=imagename+'.sc.'+str(solint)
      if os.path.exists(sc_imagename+'.mask'):
        os.system("rm -rf "+str(sc_imagename+'.mask')) 
      cleaning1_cube(vis,sc_imagename,'corrected',imagesize,cellsize,fits_dir2,scale,nterms,parallel,total_channels,freq,deltafreq)
      stat=image_stat(fits_dir2+'/'+str(sc_imagename)+'.image'+'.fits')
      beam=syn_beam(fits_dir2+'/'+str(sc_imagename)+'.image'+'.fits')
      fd.write(f"Image: {sc_imagename+'.image'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
     fd.close()
  
   if Cleaning1_sc_ap_cube==True:
     fd=open(fits_dir2+"/"+"image_statistics.txt", "a")
     cellsize='0.22arcsec'
     solints=solints_ap
     for solint in solints:
      caltable=vis+'.G.'+str(solint)+'.ap.selfcal'
      selfcal_ap(vis,caltable,solint)
      #shutil.rmtree(vis+'.G.selfcal')
      #mstransform(vis=vis,outputvis=vis[:-3]+str(solint)+'.ms',datacolumn='corrected') 
      #Vis=vis[:-3]+str(solint)+'.ms'
      sc_imagename=imagename+'.sc.ap.'+str(solint)
      if os.path.exists(sc_imagename+'.mask'):
        os.system("rm -rf "+str(sc_imagename+'.mask')) 
      cleaning1_cube(vis,sc_imagename,'corrected',imagesize,cellsize,fits_dir2,scale,nterms,parallel,total_channels,freq,deltafreq)
      stat=image_stat(fits_dir2+'/'+str(sc_imagename)+'.image'+'.fits')
      beam=syn_beam(fits_dir2+'/'+str(sc_imagename)+'.image'+'.fits')
      fd.write(f"Image: {sc_imagename+'.image'+'.fits'}\n")  # Write the image name
      fd.write(f"DR: {stat[0]}, Max: {stat[1]}, StdDev: {stat[2]}\n")  # Write the statistics below
      fd.write(f"Bmaj: {beam[0]}, bmin: {beam[1]}, bpa: {beam[2]}\n")
     fd.close()
  
  
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
