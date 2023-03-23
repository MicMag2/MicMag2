#import of MicMag2
from magnum import *
#plot module
from copy import deepcopy

#definition of simulation constants
dx,dy,dz = 2e-9,2e-9,1e-9
nx,ny,nz = 50,50,1
D = 1.35e-3
Msat = 4.3e5
#definition of the material
material = Material({'Ms':Msat,'alpha':0.5,'A':1e-11,'axis1':(0,0,1),'axis2':(0,1,0),'k_uniaxial':266676.1,'k_cubic':0})
#definition of the grid
mesh = RectangularMesh((nx,ny,nz),(dx,dy,dz),'xy')
world = World(mesh,Body("sample",material,Cuboid((0,0,0),(nx*dx,ny*dy,nz*dz))))
#Initialization of the solver with the corresponding modules, including Topology
solver = create_solver(world,[ExchangeField,StrayField,AnisotropyField,DMIField,Topology],log=True)
#DMI tensor
solver.state.Dx.fill((0, D, 0))
solver.state.Dy.fill((-D, 0, 0))
solver.state.Dz.fill((0, 0, 0))
#Initialization of the configuration
for i in range(nx):
    for j in range(ny):
        if (i-nx/2)**2+(j-ny/2)**2 <= 20:
            solver.state.M.set(i,j,0,(0,0,-Msat))
        else:
            solver.state.M.set(i,j,0,(0,0,Msat))

#Relaxation of the configuration 
solver.relax()
#Saving the configuration
writeOMF("relax.omf",solver.state.M)
#Plotting the magnetization
writeImage2("", "", r"$M$",".png", solver.state.M , "magnetisation.png", "colorhsv")

#Set the method "continuous" for the topological charge
solver.state.topo_method = "continuous"
#Output of the topological charge
print("continuous",solver.state.Q)
#Analogous approach for the other methods
solver.state.topo_method = "berg-luescher-dual-lattice"
print("berg-luescher-dual-lattice",solver.state.Q)
solver.state.topo_method = "berg-luescher"
print("berg-luescher",solver.state.Q)


#Copy of topological charge density
field = Field(mesh)
field.assign(solver.state.Q_density)
#rescaling to 1/(nm*nm) as the unit
field.scale(1e-9**2)
#Plot and save the plot for the topological charge density calculated by the Berg-Lüscher method
writeImage2("Blue-Black-Red", "", r"Q $\left[\frac{1}{\mathrm{nm}^2}\right]$",".png", field , "topologicalchargedensity.png","pixel")