import matplotlib
import matplotlib.pyplot as plt
import plot as bplt
import sys, os
import hdf5_to_dict as io
import numpy as np
font = {'size' : 16}
matplotlib.rc('font', **font)

if len(sys.argv) < 3:
  print('ERROR Format is')
  print('  snap.py [dumpfile] [variable] [coords=cart/mks] [size=40] [log=True] [vmin=-3] [vmax=3] [cmap=jet] [savefig=False]')
  sys.exit()

dfnam = sys.argv[1]
if not os.path.exists(dfnam):
  print('ERROR File ' + dfnam + ' does not exist!')
  sys.exit()

hdr = io.load_hdr(dfnam)
geom = io.load_geom(hdr)
dump = io.load_dump(dfnam)

vnam = sys.argv[2]
if not vnam in dump.keys():
  print('ERROR Variable ' + vnam + ' is not contained in dump file!')
  print('Available variables are:')
  for key in dump.keys():
    if type(dump[key]) is np.ndarray:
      if (len(dump[key].shape) == 3 and dump[key].shape[0] == hdr['N1'] and
          dump[key].shape[1] == hdr['N2'] and dump[key].shape[2] == hdr['N3']):
        print(key, end=' ')
  print('')
  sys.exit()

coords = 'cart'
for arg in sys.argv:
  if arg.startswith('coords='):
    coords = str(arg[7:])
    if coords not in ['cart', 'mks']:
      print('ERROR Coordinates ' + coords + ' are not supported!')
      print('Available coordinates are:')
      print('  cart mks')
      sys.exit()

size = 40
for arg in sys.argv:
  if arg.startswith('size='):
    size = float(arg[5:])

cmap = 'jet'
for arg in sys.argv:
  if arg.startswith('cmap='):
    cmap = str(arg[5:])

logplot = True
for arg in sys.argv:
  if arg.startswith('log='):
    if arg[4:] == 'False' or arg[4:] == 'false' or not bool(arg[4:]):
      logplot = False

savefig = False
for arg in sys.argv:
  if arg.startswith('savefig='):
    if bool(arg[8:]) == False or arg[8:] == 'false':
      savefig = False
    else:
      savefig = arg[8:]

vmin = -3; vmax = 3
if vnam == 'RHO' and logplot:
  vmin = -4; vmax = 0
for arg in sys.argv:
  if arg.startswith('vmin='):
    vmin = float(arg[5:])
  if arg.startswith('vmax='):
    vmax = float(arg[5:])

IS_3D = hdr['N3'] > 1

var = dump[vnam]
if logplot:
  var = np.log10(var)

if IS_3D:
  if coords == 'mks':
    fig, (a0, a1) = plt.subplots(1,2,gridspec_kw={'width_ratios':[1,1]}, figsize=(12,6))
  elif coords == 'cart':
    fig, (a0, a1) = plt.subplots(1,2,gridspec_kw={'width_ratios':[1,1]}, figsize=(13,7))
  ax = a0
  if coords == 'mks':
    bplt.plot_X1X2(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax, 
      cbar=False, label=vnam, ticks=None, shading='gouraud')
  elif coords == 'cart':
    bplt.plot_xz(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
      cbar=False, label=vnam, ticks=None, shading='gouraud')
    ax.set_xlim([-size,size]); ax.set_ylim([-size,size])
  ax = a1
  if coords == 'mks':
    bplt.plot_X1X3(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax, 
      cbar=True, label=vnam, ticks=None, shading='gouraud')
  elif coords == 'cart':
    bplt.plot_xy(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
      cbar=True, label=vnam, ticks=None, shading='gouraud')
    ax.set_xlim([-size,size]); ax.set_ylim([-size,size])

else:
  if coords == 'mks':
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    bplt.plot_X1X2(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
      cbar=True, label=vnam, ticks=None, shading='gouraud')
  elif coords == 'cart':
    fig, ax = plt.subplots(1, 1, figsize=(7, 10))
    bplt.plot_xz(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
      cbar=True, label=vnam, ticks=None, shading='gouraud')
    ax.set_xlim([0,size]); ax.set_ylim([-size,size])


if savefig == False:
  plt.show()
else:
  plt.savefig(savefig)

