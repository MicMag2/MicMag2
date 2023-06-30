from magnum import *
import time

mu_b = 9.27401008e-24
nx = 15
ny = 15
nz = 15
dx = 2.5e-10
dz = dy = dx
Msat = 1.85*mu_b/(dx*dy*dz)
Aex = 0.5* 5e-22/dx
alpha = 1
mu=Msat*dx*dy*dz

mat = Material({
'id': 'body',
'Ms': Msat,
'alpha': 0.5,
'J': Aex*dx,
'A': Aex,
'axis1': (0,0,1),
'axis2': (0,1,0),
'k_uniaxial': 0.0,
'k_cubic': 0.0,
'mu':mu})

mesh = RectangularMesh((nx, ny, nz), (dx, dy,dz),"xyz")
world = World(mesh, Body('body',mat,Cuboid((0,0,0),(nx*dx,ny*dx,nz*dz))))
solver = create_solver(world, [FSExternalField, TemperatureField], True, log=True, evolver = "euler", step_size=1e-15)
from scipy.constants import mu_0
solver.state.H_ext_offs = (0,0, 1.0/mu_0)

for T in range(0,12,2):
	print(T)
	solver.state.t = 0
	solver.state.kelv = Field(mesh)
	solver.state.kelv.fill(T)
	solver.state.M.fill((0,0,mu))
	solver.state.kelv_seed = time.time()
	handler = DataTableLog(f"data_{T}.odt")
	solver.addStepHandler(handler, condition.Time(1e-13))#condition.Always())
	solver.solve(condition.Time(1e-10))
	solver.removeStepHandler(handler)
