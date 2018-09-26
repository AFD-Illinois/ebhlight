################################################################################
#                                                                              #
# BASE MODULE FOR PYTHON SCRIPTING                                             # 
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True

PATHS = {}
PATHS['BASE'] = os.path.join(os.path.abspath(__file__).rsplit('/', 2)[0], '')
PATHS['CORE'] = os.path.join(PATHS['BASE'], 'core')
PATHS['SCRIPT'] = os.path.join(PATHS['BASE'], 'script')
PATHS['ANALYSIS'] = os.path.join(PATHS['BASE'], 'script', 'analysis')
PATHS['MACHINE'] = os.path.join(PATHS['BASE'], 'script', 'machine')
PATHS['PROB'] = os.getcwd()
sys.path.insert(0, PATHS['SCRIPT'])
sys.path.insert(0, PATHS['ANALYSIS'])
sys.path.insert(0, PATHS['MACHINE'])
import units
cgs = units.get_cgs()
import util
import config
import hdf5_to_dict as io

def build(PROBLEM):
  config.build(PROBLEM, PATHS)
  print(os.getcwd().split('/')[-1])

