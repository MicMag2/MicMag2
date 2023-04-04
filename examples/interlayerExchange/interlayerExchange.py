from magnum import *
from math import sqrt

tn = 15*(3e-9+4e-9)                             #thickness Pt + Ta times repeats in stack [m]
tm = 15*0.9e-9                                  #thickness Co times repeats in stack [m]
t = tm+tn                                               #cell thickness
d = 1e-9                                                #cell size [m]

Ms = 6e5*tm/t                                   #saturation magnetization [A/m]
A = 1.0e-11*tm/t                                #exchange [J/m]
D = 1.25e-3*tm/t                                #DMI [J/m**2]
Hk = 0.5                                                #in-plane saturation field [T]
K = 0.5*MU0*Ms*Ms+0.5*Ms*Hk             #K_u, perpendicular magnetic anisotropy [J/m**3]
alpha = 0.5                                             #damping constant
alpha_H = 0.35 * 15                             #Spin Hall angle times repeats
J = 1.5e11                                              #current density [A/m**2]

A_int = -1e-11                  #interlayer exchange, not enough to couple domain walls
A_int = -1e-9                   #interlayer exchange, coupled domain walls

WIRE = Material({
    'id': 'wire',
    'Ms': Ms,
    'alpha': alpha,
    'A': A,
    'axis1': (0,0,1),
    'axis2': (0,1,0),
    'k_uniaxial': K,
    'k_cubic': 0.0,
    'MST_xi': 0,
    'MST_alpha_hall': alpha_H})

mesh = RectangularMesh((200, 40, 3), (d, d, t))
(nx, ny, nz) = mesh.getNumNodes()
world = World(mesh, Body("all", WIRE))
solver = create_solver(world, [ExchangeField, StrayField, AnisotropyField, DMIField, MacroSpinTorque, InterlayerExchangeField], log=True)

for i in range(0,nx):
    for j in range(0,ny):
        for k in range(0,nz):
            if k == 0:
                if i < nx/4: solver.state.M.set(i,j,k, (0, 0, -Ms))
                else:  solver.state.M.set(i,j,k, (0, 0, Ms))
            elif k == 1:
                solver.state.M.set(i,j,k, (0, 0, 0))
                solver.state.Ms.set(i,j,k, 0)
            elif k == 2:
                if i < nx/4: solver.state.M.set(i,j,k, (0, 0, +Ms))
                else: solver.state.M.set(i,j,k, (0, 0, -Ms))

solver.state.Dx.fill((0, D, 0))
solver.state.Dy.fill((-D, 0, 0))
#solver.state.intExchPat = [(0,2,A_int)]
solver.state.intExchMatrix.fill((0,2,A_int))

solver.state.j = VectorField(mesh)
solver.state.j.fill((0,0,0))   
     
solver.solve(condition.Relaxed(100))    # VERY quick 'n' dirty. just for demonstration purposes!

print("switching current on")
solver.state.t = 0
solver.addStepHandler(OOMMFStorage("interlayer-Aint" + str(A_int*1e22) , "M"), condition.EveryNthSecond(0.5e-10))   
solver.state.j.fill((J,0,0))   
   
solver.solve(condition.Time(0.5e-9))