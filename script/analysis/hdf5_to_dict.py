import sys; sys.dont_write_bytecode = True
import units
import os
import glob
units = units.get_cgs()
SMALL = 1.e-30

supported_versions = ['bhl-release-1.0']

def h5_to_str(string):
  return string.decode()

def get_dumps_reduced(folder):
  return sorted(glob.glob(os.path.join(folder,'dump_*.h5')))

def read_scalar(dfile, name):
  if dfile[name].shape == ():
    return dfile[name][()]
  elif dfile[name].shape == (1,):
    return dfile[name][0]
  else:
    print('ERROR header value %s data format not recognized' % name)
    sys.exit()

def get_dumps_full(folder):
  import h5py
  alldumps = sorted(glob.glob(os.path.join(folder,'dump_*.h5')))
  fulldumps = []

  for fname in alldumps:
    with h5py.File(fname, 'r') as dfile:
      if 'FULL_DUMP' not in dfile:
        fulldumps.append(fname)
      elif dfile['FULL_DUMP'][0]:
        fulldumps.append(fname)

  return sorted(fulldumps)

def load_hdr(fname):
  import numpy as np
  import h5py
  dfile = h5py.File(fname, 'r')
  path = os.path.dirname(os.path.realpath(fname))  # more robust

  hdr = {}
  if 'header/version' in dfile.keys():
    hdr['VERSION'] = h5_to_str(read_scalar(dfile, 'header/version'))
  else:
    hdr['VERSION'] = None
  if hdr['VERSION'] not in supported_versions:
    print('ERROR file format ' + str(hdr['VERSION']) + ' not supported')
    sys.exit()

  hdr['PATH'] = path
  hdr['N1'] = read_scalar(dfile, '/header/n1')
  hdr['N2'] = read_scalar(dfile, '/header/n2')
  hdr['N3'] = read_scalar(dfile, '/header/n3')
  hdr['METRIC'] = h5_to_str(read_scalar(dfile, '/header/metric'))
  hdr['ELECTRONS'] = '/header/has_electrons' in dfile.keys()
  hdr['RADIATION'] = '/header/has_radiation' in dfile.keys()
  hdr['FLOORADV'] = '/header/has_flooradv' in dfile.keys()
  hdr['NVAR'] = read_scalar(dfile, '/header/n_prim')
  hdr['tf'] = read_scalar(dfile, '/header/tf')
  hdr['startx'] = np.array([0, read_scalar(dfile, '/geom/startx1'),
    read_scalar(dfile, '/geom/startx2'), read_scalar(dfile, '/geom/startx3')])
  hdr['dx'] = np.array([0, read_scalar(dfile, '/geom/dx1'), read_scalar(dfile, '/geom/dx2'),
    read_scalar(dfile, '/geom/dx3')])


  if hdr['METRIC'] == 'MKS':
    hdr['Rin'] = read_scalar(dfile, '/geom/mks/r_in')
    hdr['Rout'] = read_scalar(dfile, '/geom/mks/r_out')
    hdr['Reh'] = read_scalar(dfile, '/geom/mks/r_eh')
    hdr['Risco'] = read_scalar(dfile, '/geom/mks/r_isco')
    hdr['a'] = read_scalar(dfile, '/geom/mks/a')
    hdr['hslope'] = read_scalar(dfile, '/geom/mks/hslope')
    if hdr['RADIATION']:
      hdr['Rout_rad'] = read_scalar(dfile, '/geom/mks/r_out_rad')
  if hdr['METRIC'] == 'MMKS':
    hdr['Rin'] = read_scalar(dfile, '/geom/mmks/r_in')
    hdr['Rout'] = read_scalar(dfile, '/geom/mmks/r_out')
    hdr['Reh'] = read_scalar(dfile, '/geom/mmks/r_eh')
    hdr['Risco'] = read_scalar(dfile, '/geom/mmks/r_isco')
    hdr['a'] = read_scalar(dfile, '/geom/mmks/a')
    hdr['hslope'] = read_scalar(dfile, '/geom/mmks/hslope')
    hdr['poly_xt'] = read_scalar(dfile, '/geom/mmks/poly_xt')
    hdr['poly_alpha'] = read_scalar(dfile, '/geom/mmks/poly_alpha')
    hdr['mks_smooth'] = read_scalar(dfile, '/geom/mmks/mks_smooth')
    if hdr['RADIATION']:
      hdr['Rout_rad'] = read_scalar(dfile, '/geom/mmks/r_out_rad')

  hdr['gam'] = read_scalar(dfile, '/header/gam')
  hdr['cour'] = read_scalar(dfile, '/header/cour')

  hdr['vnams'] = [h5_to_str(vnam) for vnam in read_scalar(dfile, '/header/prim_names')]

  if hdr['ELECTRONS']:
    hdr['game'] = read_scalar(dfile, '/header/game')
    hdr['gamp'] = read_scalar(dfile, '/header/gamp')
    hdr['tptemin'] = read_scalar(dfile, '/header/tptemin')
    hdr['tptemax'] = read_scalar(dfile, '/header/tptemax')
    hdr['fel0'] = read_scalar(dfile, '/header/fel0')

  if hdr['RADIATION']:
    hdr['Mbh'] = read_scalar(dfile, '/header/Mbh')
    hdr['mbh'] = hdr['Mbh']/units['MSOLAR']
    hdr['nth'] = read_scalar(dfile, '/header/nth')
    hdr['nphi'] = read_scalar(dfile, '/header/nphi')
    hdr['maxnscatt'] = read_scalar(dfile, '/header/maxnscatt')
    hdr['numin_emiss'] = read_scalar(dfile, '/header/numin_emiss')
    hdr['numax_emiss'] = read_scalar(dfile, '/header/numax_emiss')
    hdr['nubins_emiss'] = read_scalar(dfile, '/header/nubins_emiss')
    hdr['lnumin_emiss'] = np.log(hdr['numin_emiss'])
    hdr['lnumax_emiss'] = np.log(hdr['numax_emiss'])
    hdr['dlnu_emiss'] = (hdr['lnumax_emiss'] - hdr['lnumin_emiss'])/hdr['nubins_emiss']
    hdr['nu_emiss'] = np.zeros(hdr['nubins_emiss'])
    for n in range(hdr['nubins_emiss']):
      hdr['nu_emiss'][n] = np.exp(hdr['lnumin_emiss'] + (0.5 + n)*hdr['dlnu_emiss'])
    hdr['numin_spec'] = read_scalar(dfile, '/header/numin_spec')
    hdr['numax_spec'] = read_scalar(dfile, '/header/numax_spec')
    hdr['nubins_spec'] = read_scalar(dfile, '/header/nubins_spec')
    hdr['lnumin_spec'] = np.log(hdr['numin_spec'])
    hdr['lnumax_spec'] = np.log(hdr['numax_spec'])
    hdr['dlnu_spec'] = (hdr['lnumax_spec'] - hdr['lnumin_spec'])/hdr['nubins_spec']
    hdr['nu_spec'] = np.zeros(hdr['nubins_spec'])
    for n in range(hdr['nubins_spec']):
      hdr['nu_spec'][n] = np.exp(hdr['lnumin_spec'] + (0.5 + n)*hdr['dlnu_spec'])
    hdr['thetae_max'] = read_scalar(dfile, '/header/thetae_max')
    hdr['kdotk_tol'] = read_scalar(dfile, '/header/kdotk_tol')
    hdr['L_unit'] = read_scalar(dfile, '/header/units/L_unit')
    hdr['M_unit'] = read_scalar(dfile, '/header/units/M_unit')
    hdr['T_unit'] = read_scalar(dfile, '/header/units/T_unit')
    hdr['RHO_unit'] = read_scalar(dfile, '/header/units/RHO_unit')
    hdr['U_unit'] = read_scalar(dfile, '/header/units/U_unit')
    hdr['B_unit'] = read_scalar(dfile, '/header/units/B_unit')
    hdr['Ne_unit'] = read_scalar(dfile, '/header/units/Ne_unit')
    hdr['Thetae_unit'] = read_scalar(dfile, '/header/units/Thetae_unit')
    nth = hdr['nth']
    nphi = hdr['nphi']
    hdr['dOmega'] = np.zeros([nth, nphi])
    for i in range(nth):
      for j in range(nphi):
        th0 = np.pi/nth*(i)
        th1 = np.pi/nth*(i + 1)
        hdr['dOmega'][i,j] = 2*np.pi/nphi*(np.cos(th0) - np.cos(th1))

  hdr['DTd'] = read_scalar(dfile, 'dump_cadence')
  hdr['DTf'] = read_scalar(dfile, 'full_dump_cadence')/hdr['DTd']

  if hdr['METRIC'] in ['MKS', 'MMKS']:
    hdr['Rout_vis'] = read_scalar(dfile, '/extras/r_out_vis')
    if hdr['RADIATION']:
      hdr['LEdd'] = 4.*np.pi*units['GNEWT']*hdr['Mbh']*units['MP']*units['CL']/units['THOMSON']
      hdr['nomEff'] = 0.1
      hdr['MdotEdd'] = hdr['LEdd']/(hdr['nomEff']*units['CL']**2)

  hdr['FULL_DUMP'] = read_scalar(dfile, 'is_full_dump')

  dfile.close()

  return hdr

def load_geom(hdr, recalc=True, use_3d_metrics=True, full=False):
  import numpy as np
  import h5py
  import os
  use_2d_metrics = not use_3d_metrics
  use_2d_metrics = False
  dfile = h5py.File(os.path.join(hdr['PATH'], 'grid.h5'), 'r')

  geom = {}
  
  geom['full'] = full
  
  geom['gcov'] = np.array(dfile['gcov'])
  geom['gcon'] = np.array(dfile['gcon'])
  geom['gdet'] = np.array(dfile['gdet'])
  geom['alpha'] = np.array(dfile['alpha'])

  if use_2d_metrics:
    geom['gcov'] = geom['gcov'][:,:,0]
    geom['gcon'] = geom['gcon'][:,:,0]
    geom['gdet'] = geom['gdet'][:,:,0]

  geom['X1'] = np.array(dfile['Xharm'][:,:,:,1])
  geom['X2'] = np.array(dfile['Xharm'][:,:,:,2])
  geom['X3'] = np.array(dfile['Xharm'][:,:,:,3])

  geom['x'] = np.array(dfile['Xcart'][:,:,:,1])
  geom['y'] = np.array(dfile['Xcart'][:,:,:,2])
  geom['z'] = np.array(dfile['Xcart'][:,:,:,3])

  if full:
    geom['X1f'] = np.array(dfile['XFharm'][:,:,:,1])
    geom['X2f'] = np.array(dfile['XFharm'][:,:,:,2])
    geom['X3f'] = np.array(dfile['XFharm'][:,:,:,3])

    geom['xf'] = np.array(dfile['XFcart'][:,:,:,1])
    geom['yf'] = np.array(dfile['XFcart'][:,:,:,2])
    geom['zf'] = np.array(dfile['XFcart'][:,:,:,3])

    geom['Lambda_h2cart_con'] = np.array(dfile['Lambda_h2cart_con'])
    geom['Lambda_h2cart_cov'] = np.array(dfile['Lambda_h2cart_cov'])

  if hdr['METRIC'] in ['MKS', 'MMKS']:
    geom['r'] = np.array(dfile['Xbl'][:,:,:,1])
    geom['th'] = np.array(dfile['Xbl'][:,:,:,2])
    geom['phi'] = np.array(dfile['Xbl'][:,:,:,3])
    if hdr['N3'] == 1:
      geom['phi'][:,:,:] = 0.
    geom['rcyl'] = geom['r']*np.sin(geom['th'])
    geom['rcyl'][:,0,:] = 0.
    geom['rcyl'][:,-1,:] = 0.
    if full:
      geom['Lambda_h2bl_con'] = np.array(dfile['Lambda_h2bl_con'])
      geom['Lambda_h2bl_cov'] = np.array(dfile['Lambda_h2bl_cov'])
      geom['Lambda_bl2cart_con'] = np.array(dfile['Lambda_bl2cart_con'])
      geom['Lambda_bl2cart_cov'] = np.array(dfile['Lambda_bl2cart_cov'])
      if use_2d_metrics:
        geom['Lambda_h2bl_con'] = geom['Lambda_h2bl_con'][:,:,0]
        geom['Lambda_h2bl_cov'] = geom['Lambda_h2bl_cov'][:,:,0]

  dfile.close()

  return geom

def load_dump(fname, geom=None):
  if geom == None:
    hdr = load_hdr(fname)
    geom = load_geom(hdr)

  import h5py
  import numpy as np
  hdr = load_hdr(fname)

  dfile = h5py.File(fname, 'r')
  dump = {}
  dump['hdr'] = hdr
  dump['t'] = read_scalar(dfile, 't')
  dump['dt'] = read_scalar(dfile, 'dt')
  dump['nstep'] = read_scalar(dfile, 'n_step')
  dump['ndump'] = read_scalar(dfile, 'n_dump')
  dump['is_full_dump'] = read_scalar(dfile, 'is_full_dump')

  for n in range(hdr['NVAR']):
    dump[hdr['vnams'][n]] = np.array(dfile['prims'][:,:,:,n])

  dump['jcon'] = dfile['jcon']
  if dump['is_full_dump']:
    keys = []
    keys += ['jcon', 'divb', 'fail']
    if 'Qvisc' in dump.keys():
      keys += ['Qvisc']
    if 'Qvisc_e' and 'Qvisc_p' in dump.keys():
      keys += ['Qvisc_e', 'Qvisc_p']
    if hdr['RADIATION']:
      keys += ['Nsph', 'nph', 'Nem', 'Nabs', 'Nsc', 'Jrad', 'Rmunu', 'nuLnu']
    if hdr['ELECTRONS'] and hdr['RADIATION']:
      keys += ['Qcoul']
    for key in keys:
      dump[key] = np.array(dfile[key]) + SMALL
    # divb calculator is anomalous at polar axes
    if hdr['METRIC'] == 'MKS':
      dump['divb'][:,0,:] = SMALL
      dump['divb'][:,-1,:] = SMALL
   
    if hdr['RADIATION']:
      dump['Jem'] = np.array(dfile['Jrad'][0,:,:,:]) + SMALL
      dump['Jabs'] = np.array(dfile['Jrad'][1,:,:,:]) + SMALL
      for n in range(2, hdr['maxnscatt'] + 2):
        varnam = 'Jsc%d' % (n-1)
        if n == hdr['maxnscatt'] + 1:
          varnam += '+'
        dump[varnam] = np.array(dfile['Jrad'][n,:,:,:]) + SMALL
      dump['Jsc'] = np.array(dfile['Jrad'])[2:,:,:,:].sum(axis=0) + SMALL
      dump['Jtot'] = dump['Jem'] - dump['Jabs'] + dump['Jsc']
      dump['Qem'] = dump['Nem']/hdr['DTd']*dump['UU']/dump['Jem']
      dump['Qsc'] = dump['Nsc']/hdr['DTd']*dump['UU']/dump['Jsc']
      dump['Qtot'] = (dump['Nem']+dump['Nsc'])/hdr['DTd']*dump['UU']/dump['Jtot'].clip(min=SMALL)
      if 'Nsuper' in dfile.keys():
        dump['Nsuper'] = np.array(dfile['Nsuper']) + SMALL
      if 'Esuper' in dfile.keys():
        dump['Esuper'] = np.array(dfile['Esuper']) + SMALL
  
  if hdr['ELECTRONS']:
    dump['Thetae'] = units['MP']/units['ME']*dump['KEL']*dump['RHO']**(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['up'] = dump['UU'] - dump['ue']
    dump['Thetap'] = (hdr['gamp']-1.)*dump['up']/dump['RHO']
    dump['TpTe'] = (hdr['gamp']-1.)*dump['up']/((hdr['game']-1.)*dump['ue'])
  elif hdr['RADIATION']:
    dump['Thetae'] = dump['UU']/dump['RHO']*hdr['Thetae_unit']

  dump['PRESS'] = (hdr['gam'] - 1)*dump['UU']
  dump['ENT'] = (hdr['gam']-1.)*dump['UU']*(dump['RHO']**(-1.0*hdr['gam']))
  dump['Theta'] = dump['PRESS']/dump['RHO']

  ucon, ucov, bcon, bcov = get_state(dump, geom)

  dump['ucon'] = ucon
  dump['ucov'] = ucov
  dump['bcon'] = bcon
  dump['bcov'] = bcov

  dump['bsq'] = (bcon*bcov).sum(axis=-1)

  dump['beta'] = 2.*(dump['PRESS']/(dump['bsq'] + SMALL))
  dump['sigma'] = dump['bsq']/dump['RHO']

  dump['ucon'] = ucon
  dump['ucov'] = ucov
  dump['bcon'] = bcon
  dump['bcov'] = bcov

  # WARNING jcon is in units such that it is (current density)/(4 pi)
  if 'jcon' in dfile.keys():
    jcov = np.zeros([hdr['N1'], hdr['N2'], hdr['N3'], 4])
    for mu in range(4):
      jcov[:,:,:,mu] = (dump['jcon'][:,:,:,:]*geom['gcov'][:,:,:,mu,:]).sum(axis=-1)
    dump['jcov'] = jcov
    dump['j2'] = (jcov[:,:,:,:]*dump['jcon'][:,:,:,:]).sum(axis=-1)

  if hdr['RADIATION']:
    dump['ur'] = (dump['ucon'][:,:,:,None,:]*dump['Rmunu'][:,:,:,:,:]).sum(axis=-1)
    dump['ur'] = (dump['ucov'][:,:,:,:]*dump['ur'][:,:,:,:]).sum(axis=-1)
    dump['ur'] = np.clip(dump['ur'], SMALL, None)
    dump['betar'] = dump['PRESS']/(1./3.*dump['ur'][:,:,:])

  if geom['full']:
    if hdr['METRIC'] in ['MKS', 'MMKS']:
      dump['ucon_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_con'], dump['ucon'])
      dump['ucov_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_cov'], dump['ucov'], transpose=True)
      dump['bcon_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_con'], dump['bcon'])
      dump['bcov_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_con'], dump['bcov'], transpose=True)

    dump['ucon_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'], dump['ucon'])
    dump['ucov_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_cov'], dump['ucov'], transpose=True)
    dump['bcon_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'], dump['bcon'])
    dump['bcov_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'], dump['bcov'], transpose=True)

  dfile.close()

  return dump

def load_diag(path, hdr = None):
  import numpy as np

  # load header
  if hdr is None:
    dfiles = sorted(glob.glob(os.path.join(path,'dump*.h5')))
    if len(dfiles) < 1:
      print("ERROR cannot read header. No dumps available.")
      sys.exit()
    hdr = load_hdr(dfiles[0])

  diag = {}
  dfile = np.loadtxt(os.path.join(path, 'diag.out')).transpose()
  
  # Sanity check against length-one diag files
  if not hasattr(dfile[0], "__len__"):
    print("WARNING dfile is length-one")
    return None
  
  diag['t']           = dfile[0]
  diag['rmed']        = dfile[1]
  diag['pp']          = dfile[2]
  diag['e']           = dfile[3]
  diag['adiabat rep'] = dfile[4]
  diag['u rep']       = dfile[5]
  diag['mdot']        = dfile[6]
  diag['edot']        = dfile[7]
  diag['ldot']        = dfile[8]
  diag['mass']        = dfile[9]
  diag['egas']        = dfile[10]
  diag['Phi']         = dfile[11]
  diag['phi']         = dfile[12]
  diag['jet_EM_flux'] = dfile[13]
  diag['divbmax']     = dfile[14]
  nbase = 14
  if hdr['RADIATION']:
    diag['step_made']      = dfile[nbase + 1]
    diag['step_abs']       = dfile[nbase + 2]
    diag['step_scatt']     = dfile[nbase + 3]
    diag['step_lost']      = dfile[nbase + 4]
    diag['step_rec']       = dfile[nbase + 5]
    diag['step_tot']       = dfile[nbase + 6]
    diag['step_sent']      = dfile[nbase + 7]
    diag['step_rcvd']      = dfile[nbase + 8]
    diag['step_made_all']  = dfile[nbase + 9]
    diag['step_abs_all']   = dfile[nbase + 10]
    diag['step_scatt_all'] = dfile[nbase + 11]
    diag['step_lost_all']  = dfile[nbase + 12]
    diag['step_rec_all']   = dfile[nbase + 13]
    diag['step_tot_all']   = dfile[nbase + 14]
    diag['step_sent_all']  = dfile[nbase + 15]
    diag['step_rcvd_all']  = dfile[nbase + 16]
    diag['step_fail_all']  = dfile[nbase + 17]
    diag['tune_emiss']     = dfile[nbase + 18]
    diag['tune_scatt']     = dfile[nbase + 19]
    diag['erad']           = dfile[nbase + 20]
    diag['lum']            = dfile[nbase + 21]
    diag['eff']            = dfile[nbase + 22]
    nbase += 22
    diag['Lum'] = diag['lum']*hdr['U_unit']*hdr['L_unit']**3/hdr['T_unit']
    diag['Mdot'] = diag['mdot']*hdr['M_unit']/hdr['T_unit']
    if hdr['ELECTRONS']:
      diag['num_super'] = dfile[nbase + 1]
      diag['lum_super'] = dfile[nbase + 2]
      nbase += 2
  diag['lum_eht'] = dfile[nbase + 1]
  diag['mdot_eh'] = dfile[nbase + 2]
  diag['edot_eh'] = dfile[nbase + 3]
  diag['ldot_eh'] = dfile[nbase + 4]
  diag['TIMER_UPDATE'] = dfile[nbase + 5]
  diag['TIMER_FLUXCALC'] = dfile[nbase + 6]
  diag['TIMER_FIXUP'] = dfile[nbase + 7]
  diag['TIMER_BOUND'] = dfile[nbase + 8]
  diag['TIMER_DIAG'] = dfile[nbase + 9]
  diag['TIMER_OUT'] = dfile[nbase + 10]
  diag['TIMER_MAKE'] = dfile[nbase + 11]
  diag['TIMER_PUSH'] = dfile[nbase + 12]
  diag['TIMER_INTERACT'] = dfile[nbase + 13]
  diag['TIMER_ALL'] = dfile[nbase + 14]
  if hdr['ELECTRONS']:
    diag['TIMER_ELECTRON'] = dfile[nbase + 15]
    nbase += 1

  # Ignore old data due to restarts
  ind = [0]
  t = diag['t'][0]
  for n in range(1, len(diag['t'])):
    if diag['t'][n] > t:
      t = diag['t'][n]
      ind.append(n)

  for key in diag:
    diag[key] = diag[key][ind]

  diag['hdr'] = hdr

  return diag

def get_state(dump, geom):
  import numpy as np
  hdr = dump['hdr']
  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']

  ucon = np.zeros([N1,N2,N3,4])
  ucov = np.zeros([N1,N2,N3,4])
  bcon = np.zeros([N1,N2,N3,4])
  bcov = np.zeros([N1,N2,N3,4])

  gcov = geom['gcov']
  gcon = geom['gcon']

  U1 = dump['U1']
  U2 = dump['U2']
  U3 = dump['U3']
  B1 = dump['B1']
  B2 = dump['B2']
  B3 = dump['B3']

  alpha = geom['alpha']
  #qsq = (gcov[:,:,None,1,1]*U1**2 + gcov[:,:,None,2,2]*U2**2 +
  #       gcov[:,:,None,3,3]*U3**2 + 2.*(gcov[:,:,None,1,2]*U1*U2 +
  #                                      gcov[:,:,None,1,3]*U1*U3 +
  #                                      gcov[:,:,None,2,3]*U2*U3))
  qsq = (gcov[:,:,:,1,1]*U1**2 + gcov[:,:,:,2,2]*U2**2 +
         gcov[:,:,:,3,3]*U3**2 + 2.*(gcov[:,:,:,1,2]*U1*U2 +
                                        gcov[:,:,:,1,3]*U1*U3 +
                                        gcov[:,:,:,2,3]*U2*U3))
  gamma = np.sqrt(1. + qsq)

  ucon[:,:,:,0] = gamma/alpha
  #ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,None,0,1]
  #ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,None,0,2]
  #ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,None,0,3]
  ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,:,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,:,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,:,0,3]

  for mu in range(4):
    #ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)
    ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,:,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  for mu in range(4):
    #bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)
    bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,:,mu,:]).sum(axis=-1)

  return ucon, ucov, bcon, bcov

def grid_matrix_multiply(M, v, transpose=False):
  import numpy as np
  # numpy crashes with a memory error unless
  # I use this horrible for loop. I don't know why.
  out = np.empty_like(v)
  for i in range(M.shape[0]):
    for j in range(M.shape[1]):
      for k in range(M.shape[2]):
        if len(M.shape) > 4:
          if transpose:
            out[i,j,k] = np.dot(M[i,j,k].transpose(),v[i,j,k])
          else:
            out[i,j,k] = np.dot(M[i,j,k],v[i,j,k])
        else:
          if transpose:
            out[i,j,k] = np.dot(M[i,j].transpose(),v[i,j,k])
          else:
            out[i,j,k] = np.dot(M[i,j],v[i,j,k])
  return out
