import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# GET XZ SLICE OF GRID DATA
def flatten_xz(array, hdr, flip=False):
    sign = 1.
    flat = np.zeros([2*hdr['N1'],hdr['N2']])
    for j in range(hdr['N2']):
        for i in range(hdr['N1']):
            flat[i,j] = sign*array[hdr['N1'] - 1 - i,j,hdr['N3']//2]
        for i in range(hdr['N1']):
            flat[i+hdr['N1'],j] = array[i,j,0]
    if flip:
        flat[:,0] = 0
        flat[:,-1] = 0
    return flat

# GET XY SLICE OF GRID DATA
def flatten_xy(array, hdr):
  return np.vstack((array.transpose(),array.transpose()[0])).transpose()

def plot_X1X2(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud'):
  X1 = geom['X1'][:,:,0]
  X2 = 1.-geom['X2'][:,:,0]
  mesh = ax.pcolormesh(X1, X2, var[:,:,0], cmap=cmap, vmin=vmin, vmax=vmax)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_xticklabels([]); ax.set_yticklabels([])

def plot_X1X3(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud'):
  j = dump['hdr']['N2']//2
  X1 = geom['X1'][:,j,:]
  X3 = geom['X3'][:,j,:]
  mesh = ax.pcolormesh(X1, X3, var[:,j,:], cmap=cmap, vmin=vmin, vmax=vmax)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_xticklabels([]); ax.set_yticklabels([])

def plot_xz(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True, 
  label=None, ticks=None, shading='gouraud'):
  x = geom['x']
  y = geom['y']
  z = geom['z']
  if dump['hdr']['N3'] > 1.:
    x = flatten_xz(x, dump['hdr'], flip=True)
    y = flatten_xz(y, dump['hdr'], flip=True)
    z = flatten_xz(z, dump['hdr'])
    var = flatten_xz(var, dump['hdr'])
    rcyl = np.sqrt(x**2 + y**2)
    rcyl[np.where(x<0)] *= -1
  else:
    x = x[:,:,0]
    x[:,0] = 0; x[:,-1] = 0
    y = y[:,:,0]
    z = z[:,:,0]
    var = var[:,:,0]
    rcyl = np.sqrt(x**2 + y**2)
  mesh = ax.pcolormesh(rcyl, z, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading=shading)
  circle1=plt.Circle((0,0),dump['hdr']['Reh'],color='k'); 
  ax.add_artist(circle1)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  ax.set_xlabel('x/M'); ax.set_ylabel('z/M')

def overlay_field(ax, geom, dump, NLEV=20, linestyle='-', linewidth=1,
  linecolor='k'):
  from scipy.integrate import trapz
  hdr = dump['hdr']
  N1 = hdr['N1']; N2 = hdr['N2']
  x = flatten_xz(geom['x'], hdr).transpose()
  z = flatten_xz(geom['z'], hdr).transpose()
  A_phi = np.zeros([N2, 2*N1])
  gdet = geom['gdet'][:,:,0].transpose()
  B1 = dump['B1'].mean(axis=-1).transpose()
  B2 = dump['B2'].mean(axis=-1).transpose()
  print(gdet.shape)
  for j in range(N2):
    for i in range(N1):
      A_phi[j,N1-1-i] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx'][1]) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx'][2]))
      A_phi[j,i+N1] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx'][1]) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx'][2]))
  A_phi -= (A_phi[N2//2-1,-1] + A_phi[N2//2,-1])/2.
  Apm = np.fabs(A_phi).max()
  if np.fabs(A_phi.min()) > A_phi.max():
    A_phi *= -1.
  #NLEV = 20
  levels = np.concatenate((np.linspace(-Apm,0,NLEV)[:-1], 
                           np.linspace(0,Apm,NLEV)[1:]))
  ax.contour(x, z, A_phi, levels=levels, colors=linecolor, linestyles=linestyle,
    linewidths=linewidth)

def plot_xy(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
  label=None, ticks=None, shading='gouraud'):
  hdr = dump['hdr']
  x = geom['x']
  y = geom['y']
  #x = geom['x'][:,hdr['N2']/2,:]
  #y = geom['y'][:,hdr['N2']/2,:]
  #var = var[:,dump['hdr']['N2']/2,:]
  x = flatten_xy(x[:,dump['hdr']['N2']//2,:], dump['hdr'])
  y = flatten_xy(y[:,dump['hdr']['N2']//2,:], dump['hdr'])
  var = flatten_xy(var[:,dump['hdr']['N2']//2,:], dump['hdr'])
  mesh = ax.pcolormesh(x, y, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading=shading)
  circle1=plt.Circle((0,0),dump['hdr']['Reh'],color='k'); 
  ax.add_artist(circle1)
  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, label=label, ticks=ticks)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  ax.set_xlabel('x/M'); ax.set_ylabel('y/M')
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)
 
