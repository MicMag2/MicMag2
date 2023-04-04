### SKRIPT BY THOMAS BRIAN WINKLER
### Similar to "ENERGETICS AND DYNAMICS OF A STABLE BLOCH POINT" by Winkler et al. (2023), arxiv 2303.10091
### plotly needed for plotting the 3D structure
from magnum import *
from math import *
import numpy as np
import random
import os.path
import os
import time
import sys
import re
import matplotlib.pyplot as plt
import plotly as ply
from numpy.random import default_rng
import plotly.graph_objects as go
from plotly.subplots import make_subplots

system= "BP" #chossing configurations Bloch point (BP)
nt = 50 # thickness in units of FeGe lattice constant, chose nt=50 for HB model!
nar = 100
naz1 = 30
naz2 = 20
sample ="cylindrical" # or "rectangular"]:

a=0.4679*1e-09  # lattice constant of FeGe in m
random.seed(0)
t=(naz1+naz2)*a  #total thickness in [m]
r=nar*a       #radius in [m]
dia=2*r       #diameter in [m]
d = a*(naz1+naz2)/nt  #MM mesh size in [m]
nz = nt             # number of cells in z-direction
n =int(round(dia/d)) # number of cells in x- and y-direction


mu0 = 4*np.pi*1e-7   #mu0
alpha = 0.28         #damping constant
D  =  1.58*1e-3          #DMI [J/m**2]
A  =  8.78*1e-12         #Exchange constant in [J/m]
Ms =  0.384*1e6
mu = Ms*a**3            #magnetic moment (for HB)
J  = 2*A*a
DISK = Material({
'id': 'DISK',
'Ms': Ms,
'alpha': alpha,
'A': A,
'axis1': (0,0,1),
'axis2': (0,1,0),
'k_uniaxial': 0.0,
'k_cubic': 0.0,
'l':a,
'mu':mu,
'J':J})
relax=50.0 # degree per nanosecond relax condition

mesh = RectangularMesh((n, n, nz), (d, d, d))
if sample == "rectangular":
    disk = Cuboid((0,0,0), (n*d,n*d,nz*d))
if sample == "cylindrical":
    disk = Cylinder((n*d/2, n*d/2,0.0),(n*d/2, n*d/2,nz*d), n*d/2)
world = World(mesh, Body("disk", DISK, disk))
solver = create_solver(world, module_list =[FSStrayField, FSExchangeField,  FSDMIField],
                       finescale = True,  log=True, do_precess = True) # finescale sets either MM (False) or HB mode (True).


solver.state.M.fill((0,0,mu))

def fillDMIx(field, pos):
    x,y,z = pos
    if z<(naz1*a):
        return (-D, 0, 0)
    else:
        return (D, 0, 0)

def fillDMIy(field, pos):
    x,y,z = pos
    if z<(naz1*a):
        return (0, -D, 0)
    else:
        return (0, D, 0)

def fillDMIz(field, pos):
    x,y,z = pos
    if z<(naz1*a):
        return (0,0, 0)
    else:
        return (0, 0, 0)

solver.state.Dx = fillDMIx  #fill DMI tensor
solver.state.Dy = fillDMIy  #fill DMI tensor
solver.state.Dz = fillDMIz  #fill DMI tensor
solver.state.M.normalize(mu)
solver.relax(relax)

solver.state.M.normalize(1)
numpymag = solver.state.M.to_numpy()
plt.figure(figsize=(3,3))
plt.imshow(numpymag[:,:,0,2])
plt.colorbar()
plt.title("z-magnetization bottom")
plt.savefig("HB_BP_bottomz.png")
plt.show()
plt.figure(figsize=(3,3))
plt.imshow(numpymag[:,:,0,0])
plt.colorbar()
plt.title("x-magnetization bottom")
plt.savefig("HB_BP_bottomx.png")
plt.show()
plt.figure(figsize=(3,3))
plt.imshow(numpymag[:,:,-1,2])
plt.colorbar()
plt.title("z-magnetization top")
plt.savefig("HB_BP_topz.png")
plt.show()
plt.figure(figsize=(3,3))
plt.imshow(numpymag[:,100,:,2])
plt.colorbar()
plt.title("z-magnetization slice")
plt.savefig("HB_BP_cross.png")
plt.show()

mag = solver.state.M.to_numpy()
a_nm=0.4679

x_ = np.linspace(0., 200*a_nm, int(mag.shape[0]/10))
y_ = np.linspace(0., 200*a_nm, int(mag.shape[1]/10))
z_ = np.linspace(0., 50*a_nm, int(mag.shape[2]/2))
x, y, z = np.meshgrid(x_, y_, z_, indexing='xy')
magx =mag[::10,::10,::2,0]
magy =mag[::10,::10,::2,1]
magz =mag[::10,::10,::2,2]
magzabs = np.absolute(magz)
xrav = np.ravel(x)
yrav=np.ravel(y)
zrav = np.ravel(z)
magxrav = np.ravel(magx)
magyrav = np.ravel(magy)
magzrav = np.ravel(magz)
magzravabs = np.absolute(magzrav)
deleteindicesUP = np.argwhere(magzrav <0.9) #up-magnetized
xravUP = np.delete(xrav,deleteindicesUP)
yravUP = np.delete(yrav,deleteindicesUP)
zravUP = np.delete(zrav,deleteindicesUP)
magxravUP = np.delete(magxrav,deleteindicesUP)
magyravUP = np.delete(magyrav,deleteindicesUP)
magzravUP = np.delete(magzrav,deleteindicesUP)
magzravabsUP = np.delete(magzravabs, deleteindicesUP)
deleteindicesDOWN = np.argwhere(magzrav >-0.9) #down-magnetzied
xravDOWN = np.delete(xrav,deleteindicesDOWN)
yravDOWN = np.delete(yrav,deleteindicesDOWN)
zravDOWN = np.delete(zrav,deleteindicesDOWN)
magxravDOWN = np.delete(magxrav,deleteindicesDOWN)
magyravDOWN = np.delete(magyrav,deleteindicesDOWN)
magzravDOWN = np.delete(magzrav,deleteindicesDOWN)
magzravabsDOWN = np.delete(magzravabs, deleteindicesDOWN)



trace1 = go.Isosurface(
    x=xrav,
    y=yrav,
    z=zrav,
    value=magzrav,
    isomin= -0.85,
    isomax= -0.8,
    cmin = -1,
    cmax=1,
    colorscale='RdBu',
    opacity=0.5,
    colorbar=None
        )
trace2 = go.Isosurface(
    x=xrav,
    y=yrav,
    z=zrav,
    value=magzrav,
    isomin= 0.8,
    isomax= 0.85,
    #surface_count=1,
    cmin = -1,
    cmax=1,
    colorscale='RdBu',
    opacity=0.5,
    colorbar=None
        )


trace3 = go.Isosurface(
    x=xrav,
    y=yrav,
    z=zrav,
    value=magzrav,
    isomin= -0.9,
    isomax= 0.9,
    surface_count=1,
    cmin = -1,
    cmax = 1,
    colorscale='RdBu',
    opacity=0.2,
    colorbar=None
        )
fig = make_subplots()
fig.add_trace(trace1)
fig.add_trace(trace2)
fig.add_trace(trace3)

fig.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=0.25),
                    camera_eye=dict(x=1.2, y=0.8, z=0.6),
                        xaxis_title="x (nm)",   yaxis_title= "y (nm)", zaxis_title="z (nm)",
                            zaxis= dict(tickvals=[0.0, 12.5, 25.0] )),
                          xaxis = dict(tickvals=[0, 20, 40, 60, 80]),
                      yaxis = dict(tickvals=[0, 20, 40, 60, 80])
                 )
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)
fig.update_traces(showscale=False)
fig.write_image("BP_state.png")
fig.show()
