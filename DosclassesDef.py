import toast.pipeline_tools
import numpy as np
import healpy as hp
from toast.tod import AnalyticNoise
from DosTODang import TODGb
from toast.weather import Weather
import pickle
import pandas as pd


class Dschedule:
    START_TIME = 1619830200  # saturday, 1 may 2021 0:50:00
    DURATION = 1  # 60*60*24 #(duration in seconds  )
    STOP_TIME = START_TIME + DURATION  # GMT: Wednesday, 1 January 2020 10:50:00
    SAMPLE_RATE = 1000  # Hz
    OT_ALTITUDE = 2390  # m
    OT_LON = '{}'.format(-(17 + 29 * 1 / 60 + 16 * 1 / 60 / 60))  # minus?
    OT_LAT = '{}'.format(28 + 18 * 1 / 60 + 8 * 1 / 60 / 60)
    totsamples = int((STOP_TIME - START_TIME) * SAMPLE_RATE)
    # other important variables
    ELEVATION = 60  # deg
    SCANRATE = 120  # Hz
    DET_RATE = 1000  # Hz
    weather_atacama = Weather(fname='weather_Atacama.fits')
    # weather_tenerife = Weather(fname='weather_Tenerife.fits')

#Simulation parameters
class Dsimulation:
    NSIDE = 64  # 1024
    nnz = 3  # nnz (int): the number of values per pixel.
    baseline_length = 10
    iter_max = 100
    use_noise_prior = False
    npix = 12 * NSIDE ** 2
    cmb_simu = "power_spect.txt"  # CAMBIO
    fwhm = 0.6

# File which contains distribution and characteristics of the detectors
with open('sron.pkl', 'rb') as f:
    data2 = pickle.load(f)

focalplane = data2  # this file sron.pkl already has the data calculated has 23 detectors
class Ddetector:
    detnames = list(sorted(data2.keys()))
    detquat = {x: data2[x]["quat"] for x in detnames}
    detfwhm = {x: data2[x]["fwhm_arcmin"] for x in detnames}  # 36
    detalpha = {x: data2[x]["alpha"] for x in detnames}  # 1
    detrate = {x: data2[x]["fsample"] * 10 for x in detnames}  # 1000
    detfmin = {x: data2[x]["fmin"] for x in detnames}  # 0
    detfknee = {x: data2[x]["fknee"] for x in detnames}  # 0.004
    detnet = {x: (data2[x]["NET"] / 3.2e-4) for x in detnames}  # NET 2.5
    detlabels = {x: x for x in detnames}  # separete labels 0A, 1A, 2A
    detpolcol = {x: "red" if i % 2 == 0 else "blue" for i,
    x in enumerate(detnames)}  # red or blue color to the image
def Dknee(Ddetector, fknee_select, Net_select):
    a = 0
    v = 0
    a = Ddetector.detfknee.copy()
    v = Ddetector.detnet.copy()
    for key, value in a.items():
        a[key] *= fknee_select
        v[key] *= Net_select
    return a, v


def minoise(Ddetector, fknee_select, Net_select):
    new_fknee = 0
    new_NET = 0
    new_fknee, new_NET = Dknee(Ddetector, fknee_select, Net_select)
    # create the noise using AnalyticNoise function
    noise_ar1 = AnalyticNoise(rate=Ddetector.detrate, fmin=Ddetector.detfmin, detectors=Ddetector.detnames,
                              fknee=new_fknee, alpha=Ddetector.detalpha, NET=new_NET)
    return noise_ar1


