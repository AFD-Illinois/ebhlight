import sys; sys.dont_write_bytecode=True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import util
import os
import matplotlib
import matplotlib.pyplot as plt

SADW=True

if len(sys.argv) < 3:
  print('ERROR Format is')
  print('  savg.py [/path/to/dumpfile] [variable] [SADW=True]')
  sys.exit()

path = sys.argv[1]
vnam = sys.argv[2]

hdr = io.load_hdr(path)
geom = io.load_geom(hdr)
dump = io.load_dump(path, geom)

N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']
dx1 = hdr['dx'][1]; dx2 = hdr['dx'][2]; dx3 = hdr['dx'][3]

savg = np.zeros(N1)
for i in range(N1):
  vol = 0
  vol = (dump['RHO'][i,:,:]*dx2*dx3*geom['gdet'][i,:,:]).sum(axis=-1).sum(axis=-1)
  val = (dump[vnam][i,:,:]*(dump['RHO'][i,:,:]*dx2*dx3*geom['gdet'][i,:,:])).sum(axis=-1).sum(axis=-1) 
  savg[i] = val/vol

fig, ax = plt.subplots(1,1,figsize=(12,6))
ax.plot(geom['rcyl'][:,N2//2,0], savg, color='k', linewidth=2)
ax.set_xscale('log')
ax.set_yscale('log')
plt.show()


