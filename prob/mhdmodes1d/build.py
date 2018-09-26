################################################################################
#                                                                              #
#  1D LINEAR RMHD MODES                                                        #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl
PROB = 'mhdmodes1d'

if '-idim' in sys.argv:
    IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
    IDIM = 1

assert IDIM in [1,2,3]
NTOT = 256
if IDIM == 1:
    N1 = NTOT
    N2 = 1
    N3 = 1
elif IDIM == 2:
    N1 = 1
    N2 = NTOT
    N3 = 1
elif IDIM == 3:
    N1 = 1
    N2 = 1
    N3 = NTOT

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', N1)
bhl.config.set_cparm('N2TOT', N2)
bhl.config.set_cparm('N3TOT', N3)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 5.)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('gam', 'double', default = 4./3.)
bhl.config.set_rparm('DTd', 'double', default = 5.e-1)
bhl.config.set_rparm('DTl', 'double', default = 5.e-1)
bhl.config.set_rparm('DTr', 'double', default = 10000)
bhl.config.set_rparm('DNr', 'integer', default = 100000)
bhl.config.set_rparm('nmode', 'integer', default = 2)
bhl.config.set_rparm('idim', 'integer', default = IDIM)
#bhl.config.set_rparm('cour', 'double', default=0.4)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

