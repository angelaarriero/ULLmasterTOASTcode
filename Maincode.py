import toast
import toast.pipeline_tools
from toast.mpi import MPI
import numpy as np
import healpy as hp
from toast.tod import AnalyticNoise
import os
#####-------##########
from TODangela import (TODGb)
import classesDef as cd
########----##############

#CMB INPUT MAP
MAPA_SIM= cd.cmbinput(cd.Dsimulation)
#prueba
#to change the polarization angles of the detectors
cd.Ddetector.detquat['2']=[0.027150201580442457,0.0028586627064642166,0.38116388211451296,0.9241043174734472]
cd.Ddetector.detquat['0']=[0.01545194635831245,-0.003324119502272795,0.923267816842732,0.3838316375257909]
cd.Ddetector.detquat['7']=[0.01954527152183795,-0.030073028133292462,0.922745744117246,0.3837367418602782]
cd.Ddetector.detquat['12']=[0.02913085238689356,-0.012026206551948448,0.3814051931803069,0.9238705766557758]
cd.Ddetector.detquat['10']=[0.006333093068507152,-0.0152386367540685,0.9233023087595494,0.38371932779226653]
cd.Ddetector.detquat['17']=[-0.004648629965548112,0.028130068744699158,-0.9231219361253306,-0.3834487977762909]
cd.Ddetector.detquat['22']=[-0.031104438072555272,0.026928211119763415,-0.38174173204247064,-0.9233529311131766]
cd.Ddetector.detquat['20']=[0.0027946836751364147,0.02715686177532145,-0.9232036569137106,-0.38334019155682547]

#1/F AND WHITE NOISE FUNCTION---> KNEE FRECUENCY VALUE, NET value, fmin
N2=cd.minoise(cd.Ddetector,250,1) #1/f 1Hz

weth_fil=cd.Dschedule.weather_atacama
outdir2="/scratch/aarriero/main_docs/resultados/map_maker_H_1fA_l50"
outprefix2="toast_testH_JAN3000T1fA2_l50_scan10"

Noise2,data2=cd.mifuncion(cd.Ddetector,cd.Dschedule,N2,weth_fil,cd.focalplane,
                            MAPA_SIM,cd.Dsimulation,outdir2,outprefix2)


tod_er2 = data2.obs[0]["tod"]
dets = tod_er2.local_dets[::]
sky_signal_full2=[]
full_signal2=[]
ground_sig2=[]
atmosphere_signal2=[]
signal2=[]
for det in dets:
        sky_signal_full2.append(tod_er2.local_signal(det, "sky_signal"))
        full_signal2.append(tod_er2.local_signal(det, "full_signal"))
        #atmosphere_signal2.append(tod_er2.local_signal(det, "atmosphere"))#
        signal2.append(tod_er2.local_signal(det, "signal"))
times2 = np.array(tod_er2.local_times())
detss2= np.array(dets)

#import pickle
#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/noiseH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(Noise2, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/timesH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(times2, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/skyH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(sky_signal_full2, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/fullH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(full_signal2, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/atmosphereH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(atmosphere_signal2, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/signalH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(signal2, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open('/scratch/aarriero-ext/new_simulation/Hatacama/map_maker_H_1fA_l50/dets_nameH1faz1000.pickle', 'wb') as handle:
#    pickle.dump(detss2, handle, protocol=pickle.HIGHEST_PROTOCOL)
#otro cambio

