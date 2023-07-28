from magnum import *
import magnum.magneto as magneto
import random
import os.path
import numpy as np
import os
import magnum.magneto as magneto
#from image_storage3 import ImageStorage3,writeImage3
#from image_storage2 import ImageStorage2,writeImage2
import time
from scipy.constants import mu_0

nx = 60
ny = 60
nz = 1
dx = 0.271e-9
dz = dy = dx
Msat = 1.1e6
Dind = 3.9e-3
Aex = 2e-12
Ku1 = 2.5e6
alpha = 0.05
B_ext = (0,0,1.5)
mu=Msat*dx*dy*dz

mat = Material({
'id': 'body',
'Ms': Msat,
'alpha': 0.5,
'J': Aex*dx,
'A': Aex,
'axis1': (0,0,1),
'axis2': (0,1,0),
'k_uniaxial': Ku1*dx*dy*dz,
'k_cubic': 0.0,
'mu':mu})
mesh = RectangularMesh((nx, ny, nz), (dx, dy,dz))
world = World(mesh, Body('body',mat,Cuboid((0,0,0),(nx*dx,ny*dy,nz*dz))))
solver = create_solver(world, [FSExchangeField,FSStrayField, FSAnisotropyField,FSDMIField,FSExternalField], True, log=True, do_precess=False)


q = np.array([0.01,0.01,1])
q = q/np.linalg.norm(q)
for i in range(0, nx):
    for j in range(0, ny):
      x = i*dx+0.5*dx
      y = j*dy+0.5*dy

      (cx, cy) = 0.5*nx*dx,0.5*ny*dy
      r = np.sqrt((x-cx)**2+(y-cy)**2)
      if r < 4e-9:
        solver.state.M.set(i, j, 0, tuple(-q*mu))
        #solver.state.f1.set(i, j, 0, 0)
        #solver.state.f2.set(i, j, 0, 0)
      else:
        solver.state.M.set(i, j, 0, tuple(q*mu))

solver.state.H_ext_offs = (B_ext[0]/mu_0, B_ext[1]/mu_0, B_ext[2]/mu_0) #fn = lambda t:vfield
solver.state.Dx.fill((0, -Dind, 0))
solver.state.Dy.fill((Dind, 0, 0))
solver.state.Dz.fill((0, 0, 0))

#solver.state.M = readOMF(f"relax.omf")
from scipy.constants import e

solver.solve(condition.Relaxed(0.5))

writeOMF(f"relax.omf", solver.state.M)
print(solver.state.E_stray/e)
print(solver.state.E_exch/e)
print(solver.state.E_aniso/e)
print(solver.state.E_dmi/e)
print(solver.state.E_ext/e)
print(solver.state.E_tot/e)
