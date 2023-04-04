import os
import numpy as np
import matplotlib.pyplot as plt
flist = sorted([(float((".".join(ele.split("_")[1].split(".")[:-1]))),np.genfromtxt(ele)) for ele in os.listdir("./") if ele.endswith(".odt")])
mu_b = 9.27401008e-24 
dx = 2.5e-10
dz = dy = dx
Msat = 1.85*mu_b/(dx*dy*dz)
mu=Msat*dx*dy*dz
m = np.array([np.mean([np.linalg.norm(ele[1][i,2:5])/mu for i in range(len(ele[1])//2,len(ele[1]))]) for ele in flist])
T = np.array([ele[0] for ele in flist])
plt.figure(dpi=300)
plt.plot(T,m,"ro",label="Simulation")
from scipy.constants import Boltzmann
B = 1
mu = 1.85*mu_b
Tf = np.linspace(1e-5,10,100)
x = mu*B/(Boltzmann*Tf)
#Analytical Theory from: Magnetism in Condensed Matter, Stephen Blundell, 2001
plt.plot(Tf,1/np.tanh(x)-1/x,label=r"Analytical theory $\langle m \rangle = \mathcal{L}\left(\frac{\mu B}{\mathrm{k}_B T}\right)$")
plt.xlabel(r"$T$ [K]")
plt.ylabel(r"$\langle|\vec{m}_\mathrm{average, sample}|\rangle\ [\mathrm{1}]$")
plt.title("Simulation vs. Langevin Theory of Paramagnetism")
plt.legend()
plt.savefig("plot.png")
plt.close()
