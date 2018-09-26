################################################################################
#                                                                              #
# KOMISSAROV SHOCKTUBES                                                        #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import util
import shutil

shocks = [0,1,2,3]#,4,5,6]

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'kshocks'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
util.change_cparm('N1TOT', 512, 'build.py')
call(['python', 'build.py', '-dir', TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
for shock in shocks:
  util.change_rparm('shock',str(shock),'param_template.dat')
  call(['./bhlight', '-p', 'param_template.dat'])
  shutil.move('dumps', 'dumps_%d' % shock)
os.chdir('../')

# HIGH RESOLUTION VETTED SOLUTION
soln = np.loadtxt('data/kshocks.txt')

if AUTO:
  data = {}
  data['SOLX'] = []
  data['SOLY'] = []
  data['CODEX'] = []
  data['CODEY'] = []
  for shock in shocks:
    dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps_%d' % shock + '/dump*.h5'))
    hdr = io.load_hdr(dfiles[0])
    geom = io.load_geom(hdr)
    dump = io.load_dump(dfiles[-1], geom)
    x = geom['x'][:,0,0]
    rho = dump['RHO'][:,0,0]
    u1 = dump['U1'][:,0,0]
    x_soln = soln[:,0]
    rho_soln = soln[:,2*shock+1]
    u1_soln = soln[:,2*shock+2]
    data['SOLX'].append(x_soln)
    data['SOLY'].append(rho_soln)
    data['CODEX'].append(x)
    data['CODEY'].append(rho)
  data['THRESHOLD'] = 0.05
 
  print(data.keys())
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit() 

import matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(16.18,10))

for shock in shocks:
  dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps_%d' % shock + '/dump*.h5'))
  hdr = io.load_hdr(dfiles[0])
  geom = io.load_geom(hdr)
  dump = io.load_dump(dfiles[-1], geom)
  x = geom['x'][:,0,0]
  rho = dump['RHO'][:,0,0]
  u1 = dump['U1'][:,0,0]
  x_soln = soln[:,0]
  rho_soln = soln[:,2*shock+1]
  u1_soln = soln[:,2*shock+2]

  ax = fig.add_subplot(7, 2, 2*shock + 1)
  ax.plot(x, rho, color='r')
  ax.plot(x_soln, rho_soln, color='k', linestyle='--')
  ax.set_xlim([-2,2])
  if shock < 6:
    ax.set_xticklabels([])
  ax = fig.add_subplot(7, 2, 2*shock + 2)
  ax.plot(x, u1, color='r')
  ax.plot(x_soln, u1_soln, color='k', linestyle='--')
  ax.set_xlim([-2,2])
  if shock < 6:
    ax.set_xticklabels([])

plt.savefig('kshocks.png', bbox_inches='tight')

sys.exit()

dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)

rho_code = dump['RHO'][:,0,0]
P_code = (hdr['gam']-1.)*dump['UU'][:,0,0]
u1_code = dump['U1'][:,0,0]

#rho_code = dump['RHO'][:,0,0]
#P_code   = (gam-1.)*dump['UU'][:,0,0]/(tscale*tscale)
#u_code   = dump['U1'][:,0,0]/(tscale)
#N1 = len(rho)
x_code = geom['x'][:,0,0]

if AUTO:
  data = {}
  data['SOL'] = [[x_code, rho_code]]
  data['CODE'] = [x_code, rho_code]
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit() 

# MAKE FIGURE
import matplotlib
import matplotlib.pyplot as plt
code_col = 'r'; code_ls = ''; code_mrk = '.'
#sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(4,2,1)
ax.plot(x_code, rho_code, color=code_col, linestyle=code_ls, marker=code_mrk)
#ax.plot(x, rho, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Density')
plt.xlim([-2,2])

ax = fig.add_subplot(4,2,2)
ax.plot(x_code, P_code, color=code_col, linestyle=code_ls, marker=code_mrk)
#ax.plot(x, P, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Pressure')
plt.xlim([-2,2])

ax = fig.add_subplot(4,2,3)
ax.plot(x_code, u1_code, color=code_col, linestyle=code_ls, marker=code_mrk)
#ax.plot(x, u, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Velocity')
plt.xlim([-2,2])

ax = fig.add_subplot(4,2,4)
ax.plot(x_code, P_code/rho_code, color=code_col, linestyle=code_ls, marker=code_mrk)
#ax.plot(x, P/rho, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Temperature')
plt.xlim([-2,2])

plt.subplots_adjust(wspace=0.15)

plt.savefig('kshock_%d.png' % shock, bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

