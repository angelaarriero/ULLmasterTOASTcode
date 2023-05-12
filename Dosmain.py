import numpy as np
#import healpy as hp
from toast.tod import AnalyticNoise
import os
#####-------##########
from TODangela import (TODGb)
import DosclassesDef as cd
########----##############
import pickle
#CMB INPUT MAP

Input_CMB_map= cd.cmbinput(cd.Dsimulation)

#to change the polarization angles of the detectors
cd.Ddetector.detquat['2']=[0.027150201580442457,0.0028586627064642166,0.38116388211451296,0.9241043174734472]
cd.Ddetector.detquat['0']=[0.01545194635831245,-0.003324119502272795,0.923267816842732,0.3838316375257909]
cd.Ddetector.detquat['7']=[0.01954527152183795,-0.030073028133292462,0.922745744117246,0.3837367418602782]
cd.Ddetector.detquat['12']=[0.02913085238689356,-0.012026206551948448,0.3814051931803069,0.9238705766557758]
cd.Ddetector.detquat['10']=[0.006333093068507152,-0.0152386367540685,0.9233023087595494,0.38371932779226653]
cd.Ddetector.detquat['17']=[-0.004648629965548112,0.028130068744699158,-0.9231219361253306,-0.3834487977762909]
cd.Ddetector.detquat['22']=[-0.031104438072555272,0.026928211119763415,-0.38174173204247064,-0.9233529311131766]
cd.Ddetector.detquat['20']=[0.0027946836751364147,0.02715686177532145,-0.9232036569137106,-0.38334019155682547]

#Simulated_noise---> Detector parameters, Knee freq., NET (Noise Equivalent Temperature)

PSD_Analytic_result =cd.Simulated_noise(cd.Ddetector,250,1) #1/f 1Hz

weth_fil=cd.Dschedule.weather_atacama
outdir2="/scratch/aarriero/main_docs/resultados/map_maker_test1"
outprefix2="toast_test_"

Noise1,data1=cd.TOD_generator(cd.Ddetector,cd.Dschedule,PSD_Analytic_result,weth_fil,cd.focalplane,
                            Input_CMB_map,cd.Dsimulation,outdir2,outprefix2)


#with open('noise1.pickle', 'wb') as handle:
#    pickle.dump(Noise1, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('noise1.pickle', 'rb') as handle:
#    Noise1_pkl = pickle.load(handle)