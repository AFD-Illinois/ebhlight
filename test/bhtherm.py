################################################################################
#                                                                              #
# RADIATIVE THERMAL ATMOSPHERE IN SCHWARZSCHILD SPACETIME                      #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from shutil import copyfile
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import util
import os
from os.path import join
import units
cgs = units.get_cgs()

TMP_BUILD = 'build_tmp.py'
TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'bhtherm'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

os.chdir('../prob/' + PROBLEM)

# SET RESOLUTION
copyfile('build.py', TMP_BUILD)
util.change_cparm('N1TOT', 128, TMP_BUILD) 
util.change_cparm('N2TOT', 1, TMP_BUILD) 
util.change_cparm('N3TOT', 1, TMP_BUILD) 
util.change_cparm('N1CPU', 1, TMP_BUILD)
util.change_cparm('N2CPU', 1, TMP_BUILD)
util.change_cparm('N3CPU', 1, TMP_BUILD)

# COMPILE CODE
call(['python', TMP_BUILD, '-dir', TMP_DIR])
os.remove(TMP_BUILD)
call(['cp', 'init.txt', TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
util.change_rparm('nph_per_proc', 1.e5, 'param_template.dat')
util.change_rparm('tf', 50, 'param_template.dat')
call(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)

Er = np.zeros([hdr['N1'], hdr['N2'], hdr['N3']])
Rmunu = dump['Rmunu'][:,:,:,:,:]
ucov = dump['ucov'][:,:,:,:]
ucon = dump['ucon'][:,:,:,:]
for i in range(hdr['N1']):
  for j in range(hdr['N2']):
    for k in range(hdr['N3']):
      for mu in range(4):
        for nu in range(4):
          Er[i,j,k] += ucov[i,j,k,mu]*Rmunu[i,j,k,mu,nu]*ucon[i,j,k,nu]
Er = Er.mean(axis=-1).mean(axis=-1)
Tr = (Er*hdr['U_unit']/cgs['AR'])**(1/4)
Tg = dump['UU'][:,:,:]/dump['RHO'][:,:,:]*(hdr['gam']-1)/2*cgs['MP']*cgs['CL']**2/cgs['KBOL']
Tg = Tg.mean(axis=-1).mean(axis=-1)
r_code = geom['r'][:,0,0]

soln = np.loadtxt(join('../prob', PROBLEM, 'init.txt'), skiprows=1)
r_sol = soln[:,0]
P_sol = soln[:,1]
rho_sol = soln[:,2]
T_sol = P_sol*0.5*cgs['MP']/(rho_sol*cgs['KBOL'])

if AUTO:
  data = {}
  data['SOL'] = [r_sol, T_sol]
  data['CODE'] = [r_code, Tg]
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit() 

# MAKE FIGURE
import matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(2,1,1)
ax.plot(r_code, Tr, color='r', linewidth=2, label='Radiation')
ax.plot(r_code, Tg, color='b', linewidth=2, label='Gas')
ax.plot(r_sol, T_sol, color='k', linestyle='--', linewidth=2, label='Solution')
ax.set_xlim([3.5, 20])
ax.set_ylim([1.5e12, 2.6e12])
ax.set_ylabel('<T> (K)')
ax.legend(loc=1)
ax.set_xticklabels([])
ax.grid(True, linestyle='--')

ax = fig.add_subplot(2,1,2)
T_sol_interp = np.interp(r_code, r_sol, T_sol)
ax.plot(r_code, (Tr - T_sol_interp)/T_sol_interp, color='r', linewidth=2)
ax.plot(r_code, (Tg - T_sol_interp)/T_sol_interp, color='b', linewidth=2)
ax.axhline(0, color='k', linestyle='--', linewidth=2)
ax.set_xlim([3.5, 20])
ax.set_ylim([-0.01, 0.01])
ax.set_xlabel('Rc^2/(GM)')
ax.set_ylabel('(<T> - Tsol)/Tsol')
ax.grid(True, linestyle='--')

#plt.subplots_adjust(wspace=0.15)

plt.savefig(PROBLEM + '.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

