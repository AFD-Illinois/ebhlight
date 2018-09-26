# Splitting up a folder of dump files into several large .tar.g files can be convenient for archiving

import os
import sys
import glob
import fnmatch
import tarfile

MAXSIZE = 100 # GB

if len(sys.argv) != 2:
  print('ERROR Format is')
  print('  archive_bhlight.py [folder]')
  sys.exit()

folder = os.path.join(sys.argv[1], '')
if not os.path.isdir(folder):
  print('ERROR folder ' + folder + ' is not accessible')
  sys.exit()

ddir = os.path.join(folder, 'dumps')
rdir = os.path.join(folder, 'restarts')
xdir = os.path.join(folder, 'xmf')

dumps = glob.glob(os.path.join(ddir, 'dump_*.h5'))
dump_all = glob.glob(os.path.join(ddir, '*'))
dump_other = glob.glob(os.path.join(ddir, '*.out'))
restarts_all = glob.glob(os.path.join(rdir, '*'))
xmf_all = glob.glob(os.path.join(xdir, '*'))

tarfold = folder.rsplit('/')[-2] + '_archive'
tarfold = os.path.join(os.path.dirname(os.path.dirname(folder)), tarfold)
if not os.path.exists(tarfold):
  os.makedirs(tarfold)

dsize = os.path.getsize(dumps[0])/1.e9
dumps_per_archive = int(MAXSIZE/dsize)
print('%i dumps per archive' % dumps_per_archive)

tar = tarfile.open(os.path.join(tarfold, 'restarts.tar.gz'), 'w:gz')
for name in restarts_all:
  print(name)
  tar.add(name, arcname=os.path.join('restarts', name.rsplit('/')[-1]))
tar.close()

tar = tarfile.open(os.path.join(tarfold, 'xmf.tar.gz'), 'w:gz')
for name in xmf_all:
  print(name)
  tar.add(name, arcname=os.path.join('xmf', name.rsplit('/')[-1]))
tar.close()

ndmin = 0
ndleft = len(dumps)
nd = len(dumps)
narch = 0
while ndleft > 0:
  ndmax = min(nd, ndmin+dumps_per_archive) - 1
  dumps_this_archive = dumps[ndmin:ndmax]

  if narch == 0:
    dumps_this_archive = dumps_this_archive + dump_other

  tar = tarfile.open(os.path.join(tarfold, 'dumps_%08d.tar.gz' % narch), 'w:gz')
  for name in dumps_this_archive:
    print(name)
    tar.add(name, arcname=os.path.join('dumps', name.rsplit('/')[-1]))
  tar.close()

  ndleft = ndleft - dumps_per_archive
  ndmin = ndmin + dumps_per_archive
  narch = narch + 1

print("Now do scp -r /path/to/myarchive ${ARCHIVER}:${ARCHIVE}/myarchive")

