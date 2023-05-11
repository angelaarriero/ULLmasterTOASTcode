import healpy as hp
import subprocess as sp
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
#from tqdm import tqdm


hits_map = hp.read_map("/scratch/aarriero/main_docs/resultados/map_maker_test1/toast_test_hits.fits", field=None)
hits_map[hits_map == 0] = hp.UNSEEN
destriped = hp.read_map("/scratch/aarriero/main_docs/resultados/map_maker_test1/toast_test_destriped.fits", field=None)
destriped[destriped == 0] = hp.UNSEEN
binned = hp.read_map("/scratch/aarriero/main_docs/resultados/map_maker_test1/toast_test_binned.fits", field=None)
binned[binned == 0] = hp.UNSEEN
plt.figure(figsize=[12, 12])
hp.mollview(hits_map[0], sub=[1, 3, 1], title="24h Hits Map: 1/f Noise",unit='uK')
hp.mollview(destriped[0], sub=[1, 3, 2], title="24h Destriped Map: 1/f Noise",unit='uK')
hp.mollview(binned[0], sub=[1, 3, 3], title="24h binned Map: 1/f Noise",unit='uK')
plt.savefig("/scratch/aarriero/main_docs/resultados/map_maker_test1/graficastest1.png")