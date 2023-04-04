import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann
from scipy.optimize import curve_fit
flist = sorted([(float((".".join(ele.split("_")[1].split(".")[:-1]))),np.genfromtxt(ele)) for ele in os.listdir("./") if ele.endswith(".odt")])
mu_b = 9.27401008e-24 
dx = 2.5e-10
dz = dy = dx
Msat = 2*mu_b/(dx*dy*dz)
mu=Msat*dx*dy*dz
m = np.array([np.mean([np.linalg.norm(ele[1][i,2:5])/mu for i in range(len(ele[1])//2,len(ele[1]))]) for ele in flist])
T = np.array([ele[0] for ele in flist])

fig,ax = plt.subplots(dpi=300)
J = 0.5*5e-22
ax.plot(T,m,"ro")
def f(T,Tc,beta):
    if T<Tc:
        return (1-T/Tc)**beta
    return 0.0
f = np.vectorize(f)
ot,sigma = curve_fit(f,T,m,p0=[55,1])
Tc = ot[0]
beta = ot[1]
t = np.linspace(0,100,1000)
ax.plot(t,f(t,Tc,beta),label=f"Fit 1: $T_c={Tc:.2f}\pm{sigma[0,0]**0.5:.2f}\ \mathrm{{K}}$"+"\n"+fr"$\beta={beta:.2f}\pm{sigma[1,1]**0.5:.2f}\ \mathrm{{K}}$")
#print("Tc",Tc,sigma[0,0]**0.5)
#print("Tc",beta,sigma[1,1]**0.5)
def f(T,Tc):
    if T<Tc:
        return (1-T/Tc)**0.326419
    return 0.0
f = np.vectorize(f)
ot,sigma = curve_fit(f,T,m,p0=[59])
Tc = ot[0]
#print("Tc",Tc,sigma[0,0]**0.5)
t = np.linspace(0,100,1000)
#Theoretical values from
#https://en.wikipedia.org/wiki/Ising_critical_exponents
#Christian Holm and Wolfhard Janke , Finite-size scaling study of the three-dimensional classical Heisenberg model , Physics Letters A 173 ( 1993 ) 8-12
#beta=0.6930
#1/(0.6930)*J/Boltzmann =  52.258
fig.suptitle("Curie Temperature simulation"+"\n"+r"Theory: $\beta=0.326419$, $T_c=52.26\ \mathrm{{K}}$")
ax.plot(t,f(t,Tc),label=f"Fit 1: $T_c={Tc:.2f}\pm{sigma[0,0]**0.5:.2f}\ \mathrm{{K}}$")
ax.set_xlabel(r"$T\ [K]$")
ax.set_ylabel(r"$\langle|\vec{m}_\mathrm{average, sample}|\rangle\ [\mathrm{1}]$")
ax.legend()
fig.savefig("plot.png")
plt.close(fig)