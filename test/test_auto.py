################################################################################ 
#                                                                              # 
# RUN ALL TESTS AND CHECK FOR ACCURACY                                         # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
import util
import subprocess as sp
import numpy as np
import glob
import pickle
from scipy.interpolate import interp1d as interp
import time

# DON'T RUN IF COMPUTER IS BEING USED FOR SOMETHING ELSE
import psutil
if psutil.cpu_percent() > 25.:
  sys.exit()

# EMAIL OPTIONS
TO = []
FROM = None
SUBJECT = 'BHLIGHT TESTING REPORT'
LOGNAME = 'test_auto.txt'

# SSH OPTIONS
hostname = None
port     = None
username = None

ERROR_THRESHOLD = 0.01

#os.chdir(os.path.dirname(os.path.abspath(__file__)))
#print(os.path.abspath(__file__))

# Printing
logfile = open(LOGNAME, "w")
def print2(str):
  print(str)
  logfile.write(str + '\n')

TESTS = ['bhtherm.py', 'bondi.py', 'brem.py', 'brem_mpi.py', 'kshocks.py', 
         'mhdmodes1d.py', 'sod.py', 'sod_mpi.py', 'thermalization.py',
         'thermalization_mpi.py']

SEND_REPORT = '-email' in sys.argv
FAST = '-fast' in sys.argv
VERBOSE = '-verbose' in sys.argv

print("")                                                                      
print("********************************************************************************")
print("")                                                                      
print("                                AUTOMATED TESTING")                    
print("")                                                                      
print("********************************************************************************")

DATE = time.strftime('%Y/%m/%d')
TIME = time.strftime('%H:%M:%S')
MACHINE = os.uname()[1]
HASH = 'None'
popen = sp.Popen(['git', 'show', '-s', '--format=%H'], stdout=sp.PIPE,
                   universal_newlines=True)
for line in iter(popen.stdout.readline, ""):
  HASH = line.lstrip().rstrip()
popen = sp.Popen(['git', 'branch'], stdout=sp.PIPE, universal_newlines=True)
BRANCH = 'None'
for line in iter(popen.stdout.readline, ""):
  if line[0] == '*': 
    BRANCH = line[2:].rstrip()
print2('\n  DATE:    ' + DATE)
print2('  TIME:    ' + TIME)
print2('  MACHINE: ' + MACHINE)
print2('  BRANCH:  ' + BRANCH)
print2('  COMMIT:  ' + HASH + '\n')

def name_to_args(namestring):
  """Takes a script name which may contain CLI args
  and splits out the args.

  Assumes string is generically of the form
  '<script name>.py -arg --arg'
  """
  namestring = namestring.rstrip().lstrip()
  if namestring[-3:] == '.py':
    return [namestring]
  args = namestring.split('.py ')
  args = ([args[0].lstrip().rstrip() + '.py']
          + [s.lstrip().rstrip() for s in args[1].split()])
  return args

# USE INTERPOLATION ON A (ANALYTIC SOLUTION) TO COMPARE TO B
def sanitize_array(a):
  a = np.array(a) # ensure a is a numpy array
  if len(a.shape) == 1:
    return a
  if np.prod(a.shape[1:]) > 1:
    return a
    #raise ValueError(
    #  "Array should be 1d. Array shape = {}".format(a.shape)
    #)
  return a.reshape(a.shape[0])

def L1_norm_1d(xa, ya, xb, yb):
  #xa,ya,xb,yb = [sanitize_array(a) for a in [xa,ya,xb,yb]]
  if xa[0] > xb[0]:
    xb = xb[1:]
    yb = yb[1:]
  fa = interp(xa, ya)
  norm = 0.
  nodenom = np.max(ya) <= 1e-12
  for n in range(len(xb)):
    num = np.fabs(yb[n] - fa(xb[n]))
    denom = np.fabs((yb[n] + fa(xb[n]))/2.)
    if nodenom:
      norm += num
    else:
      norm += num/denom
  return norm

def L1_norm(xa, ya, xb, yb):
  if xa is None:
    if hasattr(ya, 'shape'):
      L1s = np.zeros(len(ya))
      for n in range(len(L1s)):
        L1s[n] = 2*np.fabs(yb[n] - ya[n])/np.fabs(ya[n] + yb[n])
      return L1s.max()

  xa,ya,xb,yb = [sanitize_array(a) for a in [xa,ya,xb,yb]]
  # special case for 0d arrays
  if len(xa) == len(xb) == len(ya) == len(yb) == 1:
    if np.abs(yb[0]) <= 1e-12:
      return np.fabs(ya[0] - yb[0])
    return np.fabs((ya[0] - yb[0])/yb[0])
  
  # special case for 2d arrays, return max L1
  if xa.ndim > 1 or xb.ndim > 1 or ya.ndim > 1 or yb.ndim > 1:
    L1s = np.zeros(len(xa))
    for n in range(len(L1s)):
      L1s[n] = L1_norm_1d(xa[n], ya[n], xb[n], yb[n])/len(xb[n])
      if (np.isnan(L1s[n])):
        print(len(xb[n]))
        print(xa[n])
        print(ya[n])
        print(xb[n])
        print(yb[n])
        sys.exit()

    return L1s.max()

  return L1_norm_1d(xa, ya, xb, yb)/len(xb)

# PUT /len(xb) INTO L1_NORM_1D AND FIX THE DIVIDE BY ZERO PROBLEM IN SOME SHOCK SOLNS!!!

#  #xa,ya,xb,yb = [sanitize_array(a) for a in [xa,ya,xb,yb]]
#  if xa[0] > xb[0]:
#    xb = xb[1:]
#    yb = yb[1:]
#  fa = interp(xa, ya)
#  print('shapes...')
#  print(xa.shape)
#  print(ya.shape)
#  norm = 0.
#  nodenom = np.max(ya) <= 1e-12
#  print(len(xb))
#  for n in range(len(xb)):
#    num = np.fabs(yb[n] - fa(xb[n]))
#    print('%e ' % yb[n] + '%e ' % fa(xb[n]) + '%e' % num)
#    denom = np.fabs((yb[n] + fa(xb[n]))/2.)
#    if nodenom:
#      norm += num
#    else:
#      norm += num/denom
#  print('norm = %e' % norm)

#  return (norm/len(xb))

FAIL = False
for TEST in TESTS:
  print2('  ' + TEST[:-3])
  args = [sys.executable, TEST, '-auto']
  if FAST:
    args += ['-fast']
  popen = sp.Popen(args,
                   stdout=sp.PIPE,
                   stderr=sp.PIPE,
                   universal_newlines=True)
  for line in iter(popen.stdout.readline, ""):
    if VERBOSE:
      print2(line.rstrip())
    if line.lstrip().rstrip() == 'BUILD SUCCESSFUL':
      print2('    BUILD SUCCESSFUL')
  print2('    RUN FINISHED')
  popen.wait()

  if not os.path.isfile('data.p'):
    raise RuntimeError("Test did not succesfully complete.")
  
  with open('data.p', 'rb') as f:
    data = pickle.load(f)
    if 'SOLX' in data.keys():
      xa = data['SOLX']
      ya = data['SOLY']
      xb = data['CODEX']
      yb = data['CODEY']
    elif 'SOL' in data.keys():
      xa = data['SOL'][0]
      ya = data['SOL'][1]
      xb = data['CODE'][0]
      yb = data['CODE'][1]
    elif 'SOL_SCALAR' in data.keys():
      xa = None
      xb = None
      ya = data['SOL_SCALAR']
      yb = data['CODE_SCALAR']
      
    if 'THRESHOLD' in data.keys():
      error_threshold = data['THRESHOLD']
    else:
      error_threshold = ERROR_THRESHOLD

  norm = L1_norm(xa, ya, xb, yb)

  print('    ERROR: %.2g %%' % (100*norm))
  if norm < error_threshold:
    print('    PASS\n')
  else:
    print2('    FAIL' + '\n')
    SUBJECT += ' - FAIL'
    FAIL = True

  sp.call(['rm', 'data.p'])

logfile.close()

if SEND_REPORT:
  import paramiko
  scp = paramiko.Transport((hostname, 22))
  client = paramiko.SSHClient()
  #client.load_system_host_keys()
  client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
  k = paramiko.RSAKey.from_private_key_file('/path/to/private/ssh/key')
  client.connect(hostname, port=port, username=username, pkey=k)
  def ssh_cmd(str):
    stdin, stdout, stderr = client.exec_command(str)
  report = 'email_report.sh'
  ssh_cmd('rm ' + report)
  ssh_cmd('echo "#!/bin/bash" >> ' + report)
  ssh_cmd('echo -ne "postfix-3.3.1/bin/sendmail " >> ' + report)
  for address in TO:
    ssh_cmd('echo -ne "' + address + ' " >> ' + report)
  ssh_cmd('echo "<<EOF" >> ' + report)
  ssh_cmd('echo "subject: ' + SUBJECT + '" >> ' + report)
  ssh_cmd('echo "from: ' + FROM + '" >> ' + report)
  with open('test_auto.txt') as lfil:
    for line in lfil:
      ssh_cmd('echo "' + line.rstrip() + '" >> ' + report)
  ssh_cmd('echo "EOF" >> ' + report)
  ssh_cmd('chmod +x ' + report)
  ssh_cmd('./' + report)
  client.close()
