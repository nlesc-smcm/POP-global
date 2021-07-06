import os
import numpy
import math
from matplotlib import pyplot

numpy.random.seed(123451)

from interface import POP
from omuse.units import units,constants

from iemic_grid import depth_array,depth_levels
from bstream import barotropic_streamfunction, overturning_streamfunction,z_from_cellcenterz


if __name__=="__main__":
#prepare the plot stuff
  pyplot.ion()
  pyplot.show()

  fig, ax=pyplot.subplots(2,1,figsize=(8,8))

  dim=12

  p=POP(number_of_workers=4, channel_type="sockets", mode='test') 
  

  cwd=os.getcwd()

  d=depth_array(96,118,12)
  d=d[1:,:]

  indices=numpy.indices(d.shape)

  dmask=d  
  depth=depth_levels(dim+1, 1.8)
  dz=depth[1:]-depth[:-1]

  p.parameters.topography_option='amuse' 
  p.set_KMT(indices[0].flatten()+1,indices[1].flatten()+1, d.flat) 
  p.parameters.horiz_grid_option='internal'
  p.parameters.vert_grid_option='amuse'     
  p.parameters.vertical_layer_thicknesses=dz*(5000 | units.m)
  p.parameters.surface_heat_flux_forcing='amuse'
  p.parameters.surface_freshwater_flux_forcing='amuse'

  print(p.parameters)
  print (p.elements)
  print (p.nodes)
  input()
  '''
  #customized wind forcing

  tau_x=numpy.zeros((96,120))
  tau_y=numpy.zeros((96,120))

  forcings=p.forcings.copy()
  channel=forcings.new_channel_to(p.forcings)

  f0=8.4e-5
  beta=1.85e-11
  tau0=0.1
  
  for i in range(0,96):
    for j in range(0,120):
      cor_par=(f0 + beta*j*4e5)/f0
      tau_x[i][j]=tau0*cor_par*math.sin(2*math.pi*j/120)

  forcings.tau_x=tau_x | units.Pa
  forcings.tau_y=tau_y | units.Pa
  channel.copy_attributes(["tau_x", "tau_y"])
  '''

  print (p.elements.lat.min().in_(units.deg),p.elements.lat.max().in_(units.deg))
  print (p.elements.lon.min().in_(units.deg),p.elements.lon.max().in_(units.deg))
  input()
  print ()
  print (p.nodes.lat.min().in_(units.deg),p.nodes.lat.max().in_(units.deg))
  print (p.nodes.lon.min().in_(units.deg),p.nodes.lon.max().in_(units.deg))
  
  xmin=p.elements.lon.min().value_in(units.deg)
  xmax=p.elements.lon.max().value_in(units.deg)
  ymin=p.elements.lat.min().value_in(units.deg)
  ymax=p.elements.lat.max().value_in(units.deg)

  tnow=p.model_time
  dt=50 | units.day
  tend=tnow+(365*100 | units.day)
  t=tnow.value_in(units.day)
  t=int(t/(365))


  while tnow< tend-dt/2:
      print ("evolving to:", tnow+dt)
      p.evolve_model(tnow+dt)
      tnow=p.model_time
      t=tnow.value_in(units.day)
      t=int(t/(365))

      dlon=p.nodes[1,0].lon-p.nodes[0,0].lon
      dlat=p.nodes[0,1].lat-p.nodes[0,0].lat
      dx=2.*math.sin(dlon.value_in(units.rad))*constants.Rearth
      dy=2.*math.sin(dlat.value_in(units.rad))*constants.Rearth
  
      zc=z_from_cellcenterz(p.nodes3d[0,0,:].z)
      dz=zc[1:]-zc[:-1]
      zmin=zc.min().value_in(units.m)
      zmax=zc.max().value_in(units.m)

      psib=barotropic_streamfunction(p.nodes3d.xvel,dz,dy)
      vmax=numpy.abs(psib).max().value_in(units.Sv)

      fig, ax=pyplot.subplots(2,1,figsize=(8,8))

      im=ax[0].imshow(psib.value_in(units.Sv).T, origin="lower", cmap="seismic", \
                      vmax=vmax, vmin=-vmax, extent=[xmin,xmax,ymin,ymax],       \
                      interpolation="none")
      fig.colorbar(im,ax=ax[0],label="psib [Sv]")


      psim=overturning_streamfunction(p.nodes3d.yvel,dz,dx)
      vmax=numpy.abs(psim).max().number

      im=ax[1].imshow(psim.number.T, origin="lower", cmap="seismic", vmax=vmax,  \
                      vmin=-vmax, extent=[ymin,ymax,zmin,zmax],                  \
                      interpolation="none", aspect="auto")
      fig.colorbar(im,ax=ax[1],label="psim [Sv]")


      pyplot.savefig("psib_psib"+str(t)+".png")
      pyplot.show()
      pyplot.pause(0.5)

      pyplot.clf()

  print("done")
