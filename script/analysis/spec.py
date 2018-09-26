import matplotlib
import matplotlib.pyplot as plt
import plot as bplt
import sys
import hdf5_to_dict as io
import numpy as np
font = {'size' : 16}
matplotlib.rc('font', **font)

if len(sys.argv) < 2:
  print('ERROR Format is ')
  print('  spec.py [filename] freq=None')
  sys.exit()

fnam = sys.argv[1]
dfolder = fnam.rsplit('/', 1)[0]

freq = None
if len(sys.argv) >= 3:
  freq = float(sys.argv[2][5:])

hdr = io.load_hdr(fnam)
dump = io.load_dump(fnam)

nu = hdr['nu_spec']

if freq is not None:
  print(dump['nuLnu'].shape)
  nuLnu = dump['nuLnu'].sum(axis=0).sum(axis=-2)
  for n in range(hdr['nth']):
    print(np.interp(freq, nu, nuLnu[n,:])*4*np.pi/hdr['dOmega'][n,:].sum(axis=-1))
    print(n)

fig, ax = plt.subplots(1, 1, figsize=(10,10))
nuLnu = dump['nuLnu'].sum(axis=0).sum(axis=0).sum(axis=0)
ax.step(nu, nuLnu, where='mid', color='k')
ax.set_xscale('log'); ax.set_yscale('log')
maxL = nuLnu.max()
ax.set_ylim([1.e-10*maxL, 1.e1*maxL])
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('nuLnu (erg s^-1)')

plt.show()

