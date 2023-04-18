import healpy as hp
import subprocess as sp
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from tqdm import tqdm


destriped1 = hp.read_map("/scratch/aarriero/main_docs/resultados/map_maker_H_1fA_l50/toast_testH_JAN3000T1fA2_l50_scan10binned.fits",field=None)
destriped1[destriped1 == 0] = hp.UNSEEN

plt.figure(figsize=[12, 12])
hp.mollview(destriped1[0], sub=[3, 2, 6], title="24h Destriped Map: 1/f Noise",unit='uK')
#plt.savefig("24h_820_simu_maps.png")