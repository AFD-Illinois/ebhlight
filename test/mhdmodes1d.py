import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
from shutil import copyfile
import glob
import numpy as np
#import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis/')
import util
import hdf5_to_dict as io
TMP_DIR = 'TMP'
TMP_BUILD = 'build_tmp.py'
util.safe_remove(TMP_DIR)

AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True
if '-idim' in sys.argv:
    IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
    IDIM = 1

RES = [64, 128, 256, 512]

util.make_dir(TMP_DIR)
os.chdir('../prob/mhdmodes1d/')

copyfile('build.py', TMP_BUILD)
# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE

for n in range(len(RES)):
  util.change_cparm('N{}TOT'.format(IDIM), RES[n], TMP_BUILD)
  #util.change_cparm('RECONSTRUCTION', 'PARA', TMP_BUILD)
  call(['python', TMP_BUILD, '-dir', TMP_DIR, '-idim', str(IDIM)])
  call(['cp', os.path.join(os.getcwd(), TMP_DIR, 'bhlight'),
        '../../test/' + TMP_DIR + '/bhlight_' + str(RES[n])])
copyfile(os.path.join(os.getcwd(), TMP_DIR, 'param_template.dat'), '../../test/' +
         TMP_DIR + '/param.dat')
util.safe_remove(TMP_BUILD)
util.safe_remove(TMP_DIR)
os.chdir('../../test/')

# LOOP OVER EIGENMODES
MODES = [0, 1, 2, 3]
NAMES = ['ENTROPY', 'SLOW', 'ALFVEN', 'FAST']
NVAR = 8
VARS = ['rho', 'u', 'u1', 'u2', 'u3', 'B1', 'B2', 'B3']

amp = 1.e-4
k1 = 2.*np.pi
var0 = np.zeros(NVAR)
var0[0] = 1.
var0[1] = 1.
var0[5 + IDIM - 1] = 1.
L1 = np.zeros([len(MODES), len(RES), NVAR])
powerfits = np.zeros([len(MODES), NVAR])

for n in range(len(MODES)):
  util.change_rparm('nmode', MODES[n], TMP_DIR + '/param.dat')
  os.chdir(TMP_DIR)

  xcode = []
  ycode = []
  ysol  = []

  # EIGENMODES
  dvar = np.zeros(NVAR)
  if MODES[n] == 0: # ENTROPY
    dvar[0] = 1.
  if MODES[n] == 1: # SLOW/SOUND
    dvar[0] = 0.580429492464
    dvar[1] = 0.773905989952
    dvar[2] = -0.253320198552
  if MODES[n] == 2: # ALFVEN
    dvar[3] = 0.480384461415
    dvar[6] = 0.877058019307
  if MODES[n] == 3: # FAST
    dvar[4] = 0.480384461415
    dvar[7] = 0.877058019307
  dvar *= amp

  # permute
  if IDIM == 2:
    dvartemp = dvar[4]
    dvar[4] = dvar[3]
    dvar[3] = dvar[2]
    dvar[2] = dvartemp
    dvartemp = dvar[7]
    dvar[7] = dvar[6]
    dvar[6] = dvar[5]
    dvar[5] = dvartemp
  if IDIM == 3:
    dvartemp = dvar[2]
    dvar[2] = dvar[3]
    dvar[3] = dvar[4]
    dvar[4] = dvartemp
    dvartemp = dvar[5]
    dvar[5] = dvar[6]
    dvar[6] = dvar[7]
    dvar[7] = dvartemp

  # RUN PROBLEM FOR EACH RESOLUTION AND ANALYZE RESULT
  for m in range(len(RES)):
    #print ['./bhlight_' + str(RES[m]), '-p', 'param.dat']
    #with open(os.devnull, "w") as fnull:
    #  call(['./bhlight_' + str(RES[m]), '-p', 'param.dat'], stdout = fnull, stderr = fnull)
    call(['./bhlight_' + str(RES[m]), '-p', 'param.dat'])

    dfiles = np.sort(glob.glob('dumps/dump*.h5'))
    hdr = io.load_hdr(dfiles[-1])
    geom = io.load_geom(hdr, recalc=True)
    dump = io.load_dump(dfiles[-1], geom)
    if IDIM == 1:
      X1 = geom['x'][:,0,0]
    if IDIM == 2:
      X1 = geom['y'][0,:,0]
    if IDIM == 3:
      X1 = geom['z'][0,0,:]
    N1 = len(X1)
    dvar_code = []
    dvar_code.append(dump['RHO'].reshape(N1) - var0[0])
    dvar_code.append(dump['UU'].reshape(N1)  - var0[1])
    dvar_code.append(dump['U1'].reshape(N1)  - var0[2])
    dvar_code.append(dump['U2'].reshape(N1)  - var0[3])
    dvar_code.append(dump['U3'].reshape(N1)  - var0[4])
    dvar_code.append(dump['B1'].reshape(N1)  - var0[5])
    dvar_code.append(dump['B2'].reshape(N1)  - var0[6])
    dvar_code.append(dump['B3'].reshape(N1)  - var0[7])
    xcode.append(X1)
    ycode.append(dvar_code[0])

    files = glob.glob('dumps/*.h5')
    for f in files:
      os.remove(f)
    files = glob.glob('restarts/*.h5')
    for f in files:
      os.remove(f)
    #call(['rm', 'dumps/'])

    dvar_sol = []
    for k in range(NVAR):
      dvar_sol.append(np.real(dvar[k])*np.cos(k1*X1))
      if abs(dvar[k]) != 0.:
        L1[n][m][k] = np.mean(np.fabs(dvar_code[k] - dvar_sol[k]))
    ysol.append(dvar_sol[0])

  print('\nL1s: (MODE = %i)' % n)
  print(L1[n,:,0])
  print(L1[n,:,1])
  print(L1[n,:,2])
  print(L1[n,:,3])
  print(L1[n,:,4])
  print(L1[n,:,5])
  print(L1[n,:,6])
  print(L1[n,:,7])

  # MEASURE CONVERGENCE
  for k in range(NVAR):
    if abs(dvar[k]) != 0.:
      powerfits[n,k] = np.polyfit(np.log(RES), np.log(L1[n,:,k]), 1)[0]

  os.chdir('../')

  if not AUTO:
    # MAKE PLOTS
    fig = plt.figure(figsize=(16.18,10))

    ax = fig.add_subplot(2,1,1)
    for m in range(len(RES)):
      ax.plot(xcode[m], ycode[m], linestyle='-', marker='o', markeredgewidth=0.,
        markersize=4)
    ax.plot(xcode[-1], ysol[-1], color='k', linestyle='--')
    plt.ylim([-amp, amp])


    ax = fig.add_subplot(2,1,2)
    for k in range(NVAR):
      if abs(dvar[k]) != 0.:
        ax.plot(RES, L1[n,:,k], marker='s', label=VARS[k])

    ax.plot([RES[0]/2., RES[-1]*2.],
      10.*amp*np.asarray([RES[0]/2., RES[-1]*2.])**-2.,
      color='k', linestyle='--', label='N^-2')
    plt.xscale('log', basex=2); plt.yscale('log')
    plt.xlim([RES[0]/np.sqrt(2.), RES[-1]*np.sqrt(2.)])
    plt.xlabel('N'); plt.ylabel('L1')
    plt.title(NAMES[MODES[n]])
    plt.legend(loc=1)
    plt.savefig('mhdmodes1d_' + NAMES[MODES[n]] + '.png', bbox_inches='tight')

if AUTO:
  data = {}
  #data['CODE_SCALAR'] = [powerfits[0,0], 
  #                powerfits[1,0], powerfits[1,1], powerfits[1,2],
  #                powerfits[2,3], powerfits[2,6],
  #                powerfits[3,4], powerfits[3,7]]
  data['CODE_SCALAR'] = [powerfits[1,0], powerfits[1,1], powerfits[1,2],
                  powerfits[2,3], powerfits[2,6],
                  powerfits[3,4], powerfits[3,7]]
  data['CODE_SCALAR'] = np.array([np.mean(data['CODE_SCALAR'])])
  data['SOL_SCALAR'] = np.array([-2.]) #*np.ones(len(data['CODE_SCALAR']))#np.zeros([len(MODES), NVAR])
  data['THRESHOLD'] = 0.15
  print(data['SOL_SCALAR'])
  print(data['CODE_SCALAR'])
  import pickle
  pickle.dump(data, open('data.p', 'wb'))

# CLEAN UP
util.safe_remove(TMP_DIR)

