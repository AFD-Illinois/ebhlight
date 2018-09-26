################################################################################
#                                                                              #
#  KOMISSAROV SHOCKTUBES (Komissarov 1999 MNRAS, 303, 343)                     #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'kshocks'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 512)
bhl.config.set_cparm('N2TOT', 1)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 0.25)
bhl.config.set_rparm('dt', 'double', default = 0.25e-6)
bhl.config.set_rparm('gam', 'double', default = 4./3.)
bhl.config.set_rparm('DTd', 'double', default = 0.25e-1)
bhl.config.set_rparm('DTl', 'double', default = 0.25e-2)
bhl.config.set_rparm('DTr', 'double', default = 10.)
bhl.config.set_rparm('DNr', 'int', default = 1000000)
bhl.config.set_rparm('tscale', 'double', default = 1.e-2)
bhl.config.set_rparm('x1Min', 'double', default=-2.)
bhl.config.set_rparm('x1Max', 'double', default=2.)
bhl.config.set_rparm('shock', 'int', default=0)
bhl.config.set_rparm('cour', 'double', default=0.5)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

