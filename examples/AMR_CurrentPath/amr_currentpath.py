#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from magnum import *
import matplotlib.pyplot as plt
from PIL import Image


###helper functions to create the shape and contacts

def f_1(x):
    yf = +np.sqrt(3)*x+np.sqrt(3)/2-0.75
    return yf

def f_2(x):
    yf = -np.sqrt(3)*x+np.sqrt(3)/2-0.75
    return yf

def f_3(y):
    if y>0:
        return 1
    if y<=0:
        return 0

def draw_triangle(colors, x_mesh,y_mesh):
    for i in range(colors.shape[0]):
        for j in range(colors.shape[1]):
            x = x_mesh[i,j]
            y = y_mesh[i,j]
            f1 = f_1(x)
            f2 = f_2(x)
            f3 = f_3(y)
            if (f1<y and f2<y and f3>0):
                colors[i,j,0] = 0
                colors[i,j,1] = 0
                colors[i,j,2] = 0
            else:
                colors[i,j,0] = 255
                colors[i,j,1] = 255
                colors[i,j,2] = 255
    return colors

def draw_contacts(colors):
    for i in range(colors.shape[0]):
        for j in range(colors.shape[1]):
            if (i>186 and j<5):
                colors[i,j,0] = 255
                colors[i,j,1] = 0
                colors[i,j,2] = 0
            if (i>186 and j>194):
                colors[i,j,0] = 0
                colors[i,j,1] = 255
                colors[i,j,2] = 0
    return colors

nx=200
ny=200
x_ = np.linspace(-0.5, 0.5, nx)
y_ = np.linspace(0, 1, ny)
(x,y) = np.meshgrid(x_,y_)
colors = np.empty((nx, ny,3), dtype=np.uint8)
colors = draw_triangle(colors, x, y)
colors = draw_contacts(colors)
plt.imshow(colors)
plt.show()
im = Image.fromarray(colors)
im.save("triangle.png")
# Calculate Current Density

Msat = 4.3e5

Chrom = Material({"sigma":7.75e6}) #define custom material if necessary
n=200
mesh = RectangularMesh((n,n,1), (2.5e-9, 2.5e-9, 1e-9), periodic_bc="")
isc = ImageShapeCreator("triangle.png", mesh)
world = World(mesh, Body("contact1", Material.Au(), isc.pick("red") ), # pick pre-defined Material
              Body("contact2", Chrom, isc.pick("green") ),
              Body("triangle", Material.Py(), isc.pick("black"))
             )

solver = create_solver(world, [CurrentPath], finescale = False)


solver.U_contact = 2e-3
j = solver.state.j.to_numpy()
fig, ax = plt.subplots(1, 2, sharey = True, figsize=(6,3))
im0 = ax[0].set_title(r"$j_x$")
im0 = ax[0].imshow(j[:,:,0,0]) #row=0, col=1
plt.colorbar(im0, ax=ax[0])
im1 = ax[1].set_title(r"$j_y$")
im1 = ax[1].imshow(j[:,:,0,1]) #row=1, col=1
plt.colorbar(im1, ax=ax[1])
plt.show()

#Calculate resistance dependent on magnetization
solver.state.M = (1,0,0)
print(r"Resistance (m = m_x):", solver.state.R_contact, "Ohm")
solver.state.M = (0,1,0)
print(r"Resistance (m = m_y):", solver.state.R_contact, "Ohm")
solver.state.M = (0,0,1)
print(r"Resistance (m = m_z):",solver.state.R_contact, "Ohm")



# In[ ]:




