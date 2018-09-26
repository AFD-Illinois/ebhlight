################################################################################
#                                                                              #
#  MACHINE-SPECIFIC FUNCTIONS                                                  #
#                                                                              #
#    OPTIONS:                                                                  #
#      COMPILER   : PATH TO COMPILER EXECUTABLE                                #
#      GSL_DIR    : PATH TO GSL INSTALLATION                                   #
#      MPI_DIR    : PATH TO MPI INSTALLATION                                   #
#      HDF5_DIR   : PATH TO HDF5 INSTALLATION                                  #
#      EXECUTABLE : BINARY WRAPPER USED TO LAUNCH BHLIGHT                      #
#                                                                              #
#    MPI_DIR AND HDF5_DIR ARE NOT REQUIRED IF COMPILER HANDLES HEADERS AND     #
#    LIBRARIES FOR THESE DEPENDENCIES                                          #
#                                                                              #
################################################################################

################################################################################
# This is an example machine file that shows how you might specify settings    #
# for a machine running Debian, assuming you have installed gsl, mpi, and hdf5 #
# parallel via the package manager.                                            #
################################################################################

import util
import sys
import os

hostnames = ['debian', 'ubuntu']
flags_base = '-Wall -Werror -fdiagnostics-color -fopenmp'

def matches_host():
  host = os.uname()[1]
  for h in hostnames:
    if h in host:
      return True
  return False

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O3 -march=native'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  host['GSL_DIR']        = ''
  host['EXECUTABLE']     = 'mpirun'

  return host

