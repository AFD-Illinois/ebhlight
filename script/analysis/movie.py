#!/usr/bin/env python

import sys, os
mlocation = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(mlocation,'..'))
import util
import hdf5_to_dict as io
from snap import make_snap
from multiprocessing import Pool

# Set your ffmpeg executable here. If not available, set codec to None
codec = 'ffmpeg'
mnam_base = 'anim'

from argparse import ArgumentParser
parser = ArgumentParser(
  description='Make a movie of your simulation.')
parser.add_argument('dumpfolder',type=str,
                    help='Folder containing data')
parser.add_argument('variable',type=str,
                    help='Variable to plot')
parser.add_argument('--coords',type=str,
                    choices=['cart','mks'],default='cart',
                    help='Coordinate system. Cartesian or Modified Kerr-Schild')
parser.add_argument('-s','--size',
                    type=float,default=40,
                    help='Size of domain to plot')
parser.add_argument('-l','--log',
                    type=bool,default=True,
                    help='Log scale?')
parser.add_argument('--vmin',
                    type=float,default=-4,
                    help='Colormap lower bound')
parser.add_argument('--vmax',
                    type=float,default=0,
                    help='Colormap upper bound')
parser.add_argument('-c','--cmap',
                    type=str,default='jet',
                    help='Colormap used')

args = parser.parse_args()

dfold = util.sanitize_path(args.dumpfolder)
if not os.path.exists(dfold):
  print('ERROR Folder ' + dfnam + ' does not exist!')
  sys.exit()

tmpdir = 'FRAMES'
util.safe_remove(tmpdir)
os.mkdir(tmpdir)

dfnams = io.get_dumps_full(dfold)
hdr = io.load_hdr(dfnams[0])
geom = io.load_geom(hdr)
num_files = len(dfnams)

def make_frame(pair):
  i,d = pair
  print("frame %d/%d" % (i,num_files))
  make_snap(d,args.variable,args.coords,
            args.size,args.cmap,args.log,
            os.path.join(tmpdir,'frame_%08d.png' % i),
            args.vmin,args.vmax,
            geom=geom)

# for pair in enumerate(dfnams):
#   make_frame(pair)
p = Pool()
p.map(make_frame,enumerate(dfnams))

if codec is not None:
  from subprocess import call
  mnam = mnam_base + '_' + args.variable + '_' + args.coords + '.mp4'
  util.safe_remove(mnam)
  call([codec, '-i', os.path.join(tmpdir, 'frame_%08d.png'), mnam])
  util.safe_remove(tmpdir)
