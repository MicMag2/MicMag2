#This simulation is inspired by (no values taken over)
#https://vampire.york.ac.uk/tutorials/simulation/curie-temperature/ (accessed 24.03.22)

from magnum import *
import time

mu_b = 9.27401008e-24 
nx = 50
ny = 50
nz = 50
dx = 2.5e-10
dz = dy = dx
Msat = 2*mu_b/(dx*dy*dz)
J = 5e-22
alpha = 1
mu=Msat*dx*dy*dz

mat = Material({
'id': 'body',
'Ms': Msat,
'alpha': 0.5,
'J': J,
'axis1': (0,0,1),
'axis2': (0,1,0),
'k_uniaxial': 0.0,
'k_cubic': 0.0,
'mu':mu})

mesh = RectangularMesh((nx, ny, nz), (dx, dy,dz),"xyz")
world = World(mesh, Body('body',mat,Cuboid((0,0,0),(nx*dx,ny*dx,nz*dz))))
solver = create_solver(world, [FSExchangeField, TemperatureField], True, log=True, evolver = "euler", step_size=1e-15)

for T in range(0,100,5):
	print(T)
	solver.state.t = 0
	solver.state.kelv = Field(mesh)
	solver.state.kelv.fill(T)
	solver.state.M.fill((0,0,mu))
	solver.state.kelv_seed = time.time()
	handler = DataTableLog(f"data_{T}.odt")
	solver.addStepHandler(handler, condition.Time(1e-13))#condition.Always())
	solver.solve(condition.Time(1e-11))
	solver.removeStepHandler(handler)
