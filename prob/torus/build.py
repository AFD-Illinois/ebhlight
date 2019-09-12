import sys; sys.path.append('../../script');
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'torus'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 128)
bhl.config.set_cparm('N2TOT', 128)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MMKS')

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', True)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)
bhl.config.set_cparm('BETA_HEAT', True)
bhl.config.set_cparm('COULOMB', True)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

# RADIATION
bhl.config.set_cparm('RADIATION', False)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('NU_BINS_EMISS', 100)
bhl.config.set_cparm('NU_BINS_SPEC', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', True)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 30001)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('rin', 'double', default = 20)     # MEDIUM_DISK 10
bhl.config.set_rparm('rmax', 'double', default = 41)    # MEDIUM_DISK 20
bhl.config.set_rparm('Rout', 'double', default = 1000.)
bhl.config.set_rparm('Rout_rad', 'double', default = 40.)
bhl.config.set_rparm('gam', 'double', default = 13./9.)
bhl.config.set_rparm('DTd', 'double', default = 5.)
bhl.config.set_rparm('DTl', 'double', default = 5.e-1)
bhl.config.set_rparm('DTr', 'double', default = 1000.)
bhl.config.set_rparm('DNr', 'integer', default = 1000)
bhl.config.set_rparm('a', 'double', default = 0.9375)
bhl.config.set_rparm('mbh', 'double', default = 1.e8)
bhl.config.set_rparm('M_unit', 'double', default = 8.e23)
bhl.config.set_rparm('u_jitter', 'double', default = 0.04)
bhl.config.set_rparm('tune_emiss', 'double', 1.0)
bhl.config.set_rparm('tune_scatt', 'double', 0.1)
bhl.config.set_rparm('t0_tune_emiss', 'double', 500)
bhl.config.set_rparm('t0_tune_scatt', 'double', 500)
bhl.config.set_rparm('MAD', 'int', default = 1)
bhl.config.set_rparm('BHflux', 'double', default = 0.)
bhl.config.set_rparm('beta', 'double', default = 100.)
bhl.config.set_rparm('numin_emiss', 'double', default=1.e10) 
bhl.config.set_rparm('numax_emiss', 'double', default=1.e16) 
bhl.config.set_rparm('numin_spec', 'double', default=1.e10) 
bhl.config.set_rparm('numax_spec', 'double', default=1.e25) 
bhl.config.set_rparm('tp_over_te', 'double', default=3)
bhl.config.set_rparm('nph_per_proc', 'double', default=1.e5)
bhl.config.set_rparm('cour', 'double', default=0.9)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

