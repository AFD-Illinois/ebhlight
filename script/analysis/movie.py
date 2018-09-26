import plot as bplt
import sys, os
sys.path.append('../')
import util
import hdf5_to_dict as io
from subprocess import call

# Set your ffmpeg executable here. If not available, set codec to None
codec = 'ffmpeg'
mnam_base = 'anim'

if len(sys.argv) < 4:
  print('ERROR Format is')
  print('  movie.py [dumpfolder] [variable] [cart/mks] [size=40] [log=True] [vmin=-3] [vmax=3]')
  sys.exit()

dfold = util.sanitize_path(sys.argv[1])
if not os.path.exists(dfold):
  print('ERROR Folder ' + dfnam + ' does not exist!')
  sys.exit()

tmpdir = 'FRAMES'
util.safe_remove(tmpdir)
#util.safe_remove(mnam)
os.mkdir(tmpdir)

dfnams = io.get_dumps_full(dfold)

for n, dfnam in enumerate(dfnams):
  inst = [sys.executable, 'snap.py', dfnam]
  for arg in sys.argv[2:]:
    inst.append(arg)

  inst.append('savefig=' + os.path.join(tmpdir, 'frame_%08d.png' % n))
  print('%d' % (n+1) + ' / %d' % len(dfnams))
  call(inst)

if not codec == None:
  from subprocess import call
  mnam = mnam_base + '_' + sys.argv[2] + '_' + sys.argv[3] + '.mp4'
  util.safe_remove(mnam)
  call([codec, '-i', os.path.join(tmpdir, 'frame_%08d.png'), mnam])
  util.safe_remove(tmpdir)


