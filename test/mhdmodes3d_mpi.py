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

RES = [16, 32, 64]#, 128]

util.make_dir(TMP_DIR)
os.chdir('../prob/mhdmodes3d/')

copyfile('build_mpi.py', TMP_BUILD)
# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE

for n in xrange(len(RES)):
  util.change_cparm('N1TOT', RES[n], TMP_BUILD)
  util.change_cparm('N2TOT', RES[n], TMP_BUILD)
  util.change_cparm('N3TOT', RES[n], TMP_BUILD)
  #util.change_cparm('N1CPU', '1', TMP_BUILD)
  #util.change_cparm('N2CPU', '1', TMP_BUILD)
  #util.change_cparm('N3CPU', '1', TMP_BUILD)
  util.change_cparm('OPENMP', 'True', TMP_BUILD)
  call(['python', TMP_BUILD, '-dir', TMP_DIR])
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
MODES = [3]
NVAR = 8
VARS = ['rho', 'u', 'u1', 'u2', 'u3', 'B1', 'B2', 'B3']

amp = 1.e-4
k1 = 2.*np.pi
k2 = 2.*np.pi
k3 = 2.*np.pi
var0 = np.zeros(NVAR)
var0[0] = 1.
var0[1] = 1.
var0[5] = 1.
L1 = np.zeros([len(MODES), len(RES), NVAR])
powerfits = np.zeros([len(MODES), NVAR])

for n in xrange(len(MODES)):
  util.change_rparm('nmode', MODES[n], TMP_DIR + '/param.dat')
  os.chdir(TMP_DIR)

  # EIGENMODES
  dvar = np.zeros(NVAR)
  if MODES[n] == 0: # ENTROPY
    dvar[0] = 1.
  if MODES[n] == 1: # SLOW/SOUND
    dvar[0] = 0.556500332363
    dvar[1] = 0.742000443151
    dvar[2] = -0.282334999306
    dvar[3] = 0.0367010491491
    dvar[4] = 0.0367010491491
    dvar[5] = -0.195509141461
    dvar[6] = 0.0977545707307
    dvar[7] = 0.0977545707307
  if MODES[n] == 2: # ALFVEN
    dvar[3] = -0.339683110243
    dvar[4] = 0.339683110243
    dvar[6] = 0.620173672946
    dvar[7] = -0.620173672946
  if MODES[n] == 3: # FAST
    dvar[0] = 0.481846076323
    dvar[1] = 0.642461435098
    dvar[2] = -0.0832240462505
    dvar[3] = -0.224080007379
    dvar[4] = -0.224080007379
    dvar[5] = 0.406380545676
    dvar[6] = -0.203190272838
    dvar[7] = -0.203190272838
  dvar *= amp
  
  # RUN PROBLEM FOR EACH RESOLUTION AND ANALYZE RESULT
  for m in xrange(len(RES)):
    #print ['./bhlight_' + str(RES[m]), '-p', 'param.dat']
    NPROC = 8
    import psutil
    NCORE = psutil.cpu_count(logical=False)
    os.environ['OMP_NUM_THREADS'] = '%d' % (NCORE/NPROC)
    call(['mpirun', '-np', '8', './bhlight_' + str(RES[m]), '-p', 'param.dat'])

    dfiles = np.sort(glob.glob('dumps/dump*.h5'))
    hdr = io.load_hdr(dfiles[-1])
    geom = io.load_geom(hdr, recalc=True)
    dump = io.load_dump(dfiles[-1], geom) 
    X1 = geom['x'][:,:,:]
    X2 = geom['y'][:,:,:]
    X3 = geom['z'][:,:,:]
    dvar_code = []
    dvar_code.append(dump['RHO'][:,:,:] - var0[0]) 
    dvar_code.append(dump['UU'][:,:,:]  - var0[1])
    dvar_code.append(dump['U1'][:,:,:]  - var0[2])
    dvar_code.append(dump['U2'][:,:,:]  - var0[3])
    dvar_code.append(dump['U3'][:,:,:]  - var0[4])
    dvar_code.append(dump['B1'][:,:,:]  - var0[5])
    dvar_code.append(dump['B2'][:,:,:]  - var0[6])
    dvar_code.append(dump['B3'][:,:,:]  - var0[7])

    #dvar_sol = []
    dvar_sol = np.zeros([RES[m], RES[m], RES[m]])
    for k in xrange(NVAR):
      #dvar_sol.append(np.real(dvar[k])*np.cos(k1*X1 + k2*X2))
      if abs(dvar[k]) != 0.:
        for i in xrange(RES[m]):
          for j in xrange(RES[m]):
            for kk in xrange(RES[m]):
              dvar_sol[i,j,kk] = np.real(dvar[k])*np.cos(k1*X1[i,j,kk] + 
                                                         k2*X2[i,j,kk] + 
                                                         k3*X3[i,j,kk])
              L1[n][m][k] = np.mean(np.fabs(dvar_code[k][i,j,kk] - 
                                            dvar_sol[i,j,kk]))

  # MEASURE CONVERGENCE
  for k in xrange(NVAR):
    if abs(dvar[k]) != 0.:
      powerfits[n,k] = np.polyfit(np.log(RES), np.log(L1[n,:,k]), 1)[0]
  
  os.chdir('../')

  if not AUTO:
    # MAKE PLOTS
    fig = plt.figure(figsize=(16.18,10))

    ax = fig.add_subplot(1,1,1)
    for k in xrange(NVAR):
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
    plt.savefig('mhdmodes3d_mpi_' + NAMES[MODES[n]] + '.png', bbox_inches='tight')

if AUTO:
  data = {}
  data['SOL'] = -2.*np.zeros([len(MODES), NVAR])  
  data['CODE'] = powerfits
  import pickle
  pickle.dump(data, open('data.p', 'wb'))

# CLEAN UP
util.safe_remove(TMP_DIR)

