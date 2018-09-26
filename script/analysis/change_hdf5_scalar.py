import h5py
import sys
from os.path import isfile
import numpy as np

if len(sys.argv) != 4:
  print('ERROR format is')
  print('  change_scalar.py [filename] [dataset name] [new value]')
  sys.exit()

fnam = sys.argv[1]
if not isfile(fnam):
  print('ERROR cannot open file %s' % fnam)
  sys.exit()

hfil = h5py.File(fnam, 'r+')
dsetnam = sys.argv[2]
dset_exists = dsetnam in hfil
if not dset_exists:
  print('ERROR cannot open dataset %s' % dsetnam)
  print('  Datasets in file:')
  #print(list(hfil.keys()))
  def print_attrs(name, obj):
    print('   ', name)
  hfil.visititems(print_attrs)
  hfil.close()
  sys.exit()

if len(hfil[dsetnam].shape) is not 1 or hfil[dsetnam].shape[0] is not 1:
  print('ERROR cannot only modify size-1 1D datasets')
  hfile.close()
  sys.exit()

print('Dataset', hfil[dsetnam])
print('Old value', hfil[dsetnam][0])
dtype = hfil[dsetnam].dtype

# Change dataset to input value
if 'float' in dtype.name:
  val = float(sys.argv[3])
elif 'int' in dtype.name:
  val = int(sys.argv[3])
else:
  print('ERROR datatype', dtype, 'not supported')
  hfil.close()
  sys.exit()
hfil[dsetnam][0] = val
print('New value', hfil[dsetnam][0])

hfil.close()

