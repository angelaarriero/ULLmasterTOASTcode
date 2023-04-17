import toast.pipeline_tools
import numpy as np
import healpy as hp
from toast.tod import AnalyticNoise
from TODangela import TODGb
from toast.weather import Weather
import pickle
import pandas as pd

class Dschedule:
    START_TIME = 1619830200 #saturday, 1 may 2021 0:50:00
    DURATION = 1 # 60*60*24 #(duration in seconds  )
    STOP_TIME = START_TIME + DURATION #GMT: Wednesday, 1 January 2020 10:50:00
    SAMPLE_RATE = 1000 #Hz
    OT_ALTITUDE = 2390 # m
    OT_LON = '{}'.format(-(17 + 29*1/60 + 16*1/60/60)) #minus?
    OT_LAT = '{}'.format(28 + 18*1/60 + 8*1/60/60)
    totsamples = int((STOP_TIME - START_TIME)*SAMPLE_RATE)
     #other important variables
    ELEVATION =60 #deg
    SCANRATE =120 #Hz
    DET_RATE =1000 #Hz
    weather_atacama = Weather(fname='weather_Atacama.fits')
    #weather_tenerife = Weather(fname='weather_Tenerife.fits')
#File which contains distribution and characteristics of the detectors  
with open('sron.pkl','rb') as f:
    data2=pickle.load(f)

focalplane=data2 #this file sron.pkl already has the data calculated has 23 detectors

class Ddetector:
    detnames = list(sorted(data2.keys()))
    detquat = {x: data2[x]["quat"] for x in detnames}
    detfwhm = {x: data2[x]["fwhm_arcmin"] for x in detnames}#36
    detalpha = {x: data2[x]["alpha"] for x in detnames}#1
    detrate = {x: data2[x]["fsample"]*10 for x in detnames}#1000
    detfmin = {x: data2[x]["fmin"] for x in detnames} #0
    detfknee = {x: data2[x]["fknee"] for x in detnames}#0.004
    detnet = {x: (data2[x]["NET"]/3.2e-4) for x in detnames} # NET 2.5
    detlabels = {x: x for x in detnames} #separete labels 0A, 1A, 2A
    detpolcol = {x: "red" if i % 2 == 0 else "blue" for i, 
                 x in enumerate(detnames)} #red or blue color to the image   
class Dsimulation:
    NSIDE=64 #1024
    nnz=3 #nnz (int): the number of values per pixel.
    baseline_length=10
    iter_max=100
    use_noise_prior=False
    npix = 12 * NSIDE ** 2
    cmb_simu="power_spect.txt" #CAMBIO
    fwhm=0.6

def Dknee(Ddetector,number,n2):
    a=0
    v=0
    a=Ddetector.detfknee.copy()
    v=Ddetector.detnet.copy()
    for key, value in a.items():
        a[key]*=number
        v[key]*=n2
    return a,v

def minoise(Ddetector,numberD,n2):
    b=0
    c=0
    b,c= Dknee(Ddetector,numberD,n2)
      #create the noise using AnalyticNoise function
    noise_ar1=AnalyticNoise(rate=Ddetector.detrate ,fmin=Ddetector.detfmin, detectors=Ddetector.detnames,
                       fknee=b , alpha=Ddetector.detalpha ,NET=c)
    return noise_ar1  

def mifuncion(Ddetector,Dschedule,noise_ar,weather_atacama,focalplane,MAPA_SIM,Dsimulation,               outdirD,outprefixD):
  #communication MPI
  data=0
  mpiworld, procs, rank = toast.mpi.get_world()
  comm = toast.mpi.Comm(mpiworld)
  obs = {}
  noisew=0
  det1=0
  obs["noise"] = noise_ar
  obs['weather']= weather_atacama
  obs['focalplane']= focalplane
  obs["site_id"] = 123
  obs["telescope_id"] = 1234
  obs["altitude"] = Dschedule.OT_ALTITUDE
  obs["fpradius"] = 3 #radius = 1.1558552389176189 deg
  obs["id"] = int(Dschedule.START_TIME * 10000)
  obs["tod"] = todgb = TODGb(comm.comm_group,
                              Ddetector.detquat,
                              Dschedule.totsamples,
                              firsttime=Dschedule.START_TIME,
                              el=Dschedule.ELEVATION,
                              site_lon=Dschedule.OT_LON,
                              site_lat=Dschedule.OT_LAT,
                              site_alt=Dschedule.OT_ALTITUDE,
                              scanrate=Dschedule.SCANRATE,
                              rate=Dschedule.DET_RATE)
  data = toast.Data(comm) 
  data.obs.append(obs)
  #name = "signal_D"
  toast.tod.OpCacheClear("signal").exec(data)
  toast.todmap.OpPointingHpix(nside=Dsimulation.NSIDE, nest=True, mode="IQU").exec(data)
  
  distmap = toast.map.DistPixels(
    data,
    nnz=Dsimulation.nnz,
    dtype=np.float32,
    )
  distmap.read_healpix_fits(MAPA_SIM)
  toast.todmap.OpSimScan(input_map=distmap, out="signal").exec(data)
  # Copy the sky signal
  toast.tod.OpCacheCopy(input="signal", output="sky_signal", force=True).exec(data)
  # Simulate noise
  toast.tod.OpSimNoise(out="signal", realization=0).exec(data)
  toast.tod.OpCacheCopy(input="signal", output="full_signal", force=True).exec(data)

  mapmaker = toast.todmap.OpMapMaker(
    nside=Dsimulation.NSIDE,
    nnz=Dsimulation.nnz,
    name="signal",
    outdir=outdirD,
    outprefix=outprefixD,
    baseline_length=Dsimulation.baseline_length,
    iter_max=Dsimulation.iter_max,
    use_noise_prior=Dsimulation.use_noise_prior,

    )
  mapmaker.exec(data)
  return noise_ar,data

def cmbinput(Dsimulation):
    input_cl = pd.read_csv(Dsimulation.cmb_simu,
                       delim_whitespace=True, index_col=0)
    lmax = input_cl.index[-1]
    cl = input_cl.divide(input_cl.index * (input_cl.index+1) / (np.pi*2), axis="index")
    #cl /= 1e12
    cl = cl.reindex(np.arange(0, lmax+1))
    cl = cl.fillna(0)
    seed = 58
    np.random.seed(seed)
    alm = hp.synalm((cl.TT, cl.EE, cl.BB, cl.TE), lmax=lmax, new=True)
    high_nside = Dsimulation.NSIDE
    cmb_map = hp.alm2map(alm, nside=high_nside, lmax=lmax,fwhm=np.radians(Dsimulation.fwhm))
    #hp.mollview(cmb_map[0], min=-300*1e-6, max=300*1e-6, unit="K", title="CMB Temperature")
    hp.mollview(cmb_map[0], title="CMB Temperature",unit='uK')
    hp.write_map("/scratch/aarriero/main_docs/ULLmasterTOASTcode/sim_map_d.fits", hp.reorder(cmb_map, r2n=True), nest=True, overwrite=True)
    MAPA_SIM="sim_map_d.fits"
    return MAPA_SIM

#prueba borrar
