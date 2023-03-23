from __future__ import division
from magnum.solver import StepHandler
import math
# import dill as pickle
import time


import magnum.evolver as evolver
#import magnum.fEvolver as fEvolver
import magnum.solver.condition as condition
import magnum.solver as solver
import magnum.magneto as magneto
import magnum.module as module
from magnum.mesh import RectangularMesh
from magnum.micromagnetics.world import *
# from constants import PI
from magnum.create_solver import create_solver
from magnum.heisenberg_model import *
from magnum.mesh import *
import magnum.solver.condition as condition
import numpy as np
import skimage.transform as strafo
import scipy as sc
import scipy as sc
from magnum.micromagnetics import StrayFieldCalculator
import multiprocessing
#import skimage.transform as strafo
import threading
#from magnum import *

# from numba import jit

class iterable_CC(type):

    def __iter__(cls):
        return iter(cls.self)


class ContinuityCheck(StepHandler):
    __metaclass__ = iterable_CC

    def __init__(self, coarseSolver, fineSolver, r_0, r_1, FinenodesPerCoarseCell, folder, Tracking, TrackingCondition,
                 fineMagIsGiven, DMI=False, mintrack=0.05, multiscaleLog=True, growshrink = False,  catchCoarseField = True, useMMstrayField = False):

        super(ContinuityCheck, self).__init__()
        self.useMMstrayField =useMMstrayField
        self.DMI = DMI
        self.solver = coarseSolver
        self.world = coarseSolver.world
        self.mesh = coarseSolver.world.mesh
        self.NumNodes = self.mesh.getNumNodes()
        self.Delta = self.mesh.getDelta()
        self.mintrack = mintrack
        self.catchCoarseField = catchCoarseField
        self.growshrink = growshrink
        self.fineSolver = fineSolver
        self.fineWorld = fineSolver.world
        self.fineMesh = fineSolver.world.mesh
        self.fineNodes = self.fineMesh.getNumNodes()
        self.l = self.fineMesh.getDelta()
        self.J = fineSolver.world.bodies[0].material.J
        # print(self.fineMesh.getDelta())
        # made for the skyrmion edge trackingte

        # if(hasattr(self, 'TrackingCond')):
        #print("self.mesh"+str(self.mesh))
        #        if not hasattr()


        self.r_0 = r_0
        self.r_1 = r_1
        self.body = self.world.bodies[0]
        self.Ms = self.body.material.Ms
        self.mu = self.body.material.mu
        if self.Ms == 0 or self.mu == 0:
            self.body = self.world.bodies[1]
            self.Ms = self.body.material.Ms
            self.mu = self.body.material.mu
        self.FinenodesPerCoarseCell = FinenodesPerCoarseCell

        self.stray = StrayFieldCalculator(self.mesh)
        xnum = int(self.fineNodes[0]/self.FinenodesPerCoarseCell[0])
        ynum = int(self.fineNodes[1]/self.FinenodesPerCoarseCell[1])
        znum = int(self.fineNodes[2]/self.FinenodesPerCoarseCell[2])
        self.stray2 = StrayFieldCalculator(RectangularMesh((xnum, ynum, znum), self.Delta))


        self.zoomsize = [r_1[0] - r_0[0], r_1[1] - r_0[1], r_1[2] - r_0[2]]
        self.oldcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.newcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.t = self.solver.state.t

        self.handletime = 0
        self.straytime = 0
        self.exchangetime = 0
        self.solvetime = 0

        self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
        self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
        self.Coarseshift = [0, 0, 0]

        if self.NumNodes[2] == 1:
            if self.NumNodes[1] == 1:
                self.dim = 1
        if self.NumNodes[2] == 1:
            if self.NumNodes[1] == 1:
                self.dim = 1
            else:
                self.dim = 2
        else:
            self.dim = 3
        self.Tracking = Tracking
        self.TrackingCond = TrackingCondition
        self.Fineshift = [0, 0, 0]
        self.fineMagIsGiven = fineMagIsGiven
        self.multiscaleLog = multiscaleLog
        self.folder = folder


    def reInit(self,tmpr0, tmpr1, tmpfineSolver, mancopy, tmpbuffer1=0, tmpbuffer2=0):
        for sh, c in self.fineSolver._Solver__step_handlers:
            tmpfineSolver._Solver__step_handlers.append((sh, c))
            self.fineSolver.removeStepHandler(sh)

        print("old field value before copy: "+str(self.fineSolver.state.fineExternalField.get(0,0,0)))
        print("tmpfineSolver field value before copy: "+str(tmpfineSolver.state.fineExternalField.get(0,0,0)))
        print("WARNING: GROW/SHRINK THE FINESOLVER ONLY (AND ONLY) WORKS FOR 1 SINGLE UNIFORM FIELD IN FineExternalField")
        self.r_0 = tmpr0
        self.r_1 = tmpr1
        self.fineWorld = tmpfineSolver.world
        self.fineMesh = tmpfineSolver.world.mesh
        self.fineNodes = tmpfineSolver.mesh.getNumNodes()
        self.l = self.fineMesh.getDelta()
        self.body = self.world.bodies[0]
        #tmpfineSolver.stop_condition.get_time_of_interest(tmpfineSolver.__state) = self.fineSolver.stop_condition.get_time_of_interest(self.fineSolver.__state)
        tmpfineSolver.state.fineExternalField.fill(mancopy)
        self.fineSolver = tmpfineSolver

        self.fineSolver.state.fineExternalField.fill(mancopy)
        print("new fineSolver field value: "+str(self.fineSolver.state.fineExternalField.get(0,0,0)))

        self.zoomsize = [self.r_1[0] - self.r_0[0], self.r_1[1] - self.r_0[1], self.r_1[2] - self.r_0[2]]
        self.oldcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.newcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.t = self.solver.state.t
        # self.fineSolver.h_try = self.solver.h_try
        self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
        self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
        # self.Coarseshift=[0,0,0]
        self.fineSolver = tmpfineSolver
        self.catchCoarseField = True

    def oneDcenter(self, A, idx):
        return (A - self.r_0[idx]) * self.FinenodesPerCoarseCell[idx] + int(self.FinenodesPerCoarseCell[idx] / 2)

    def center(self, A):
        return [self.oneDcenter(A[0], 0), self.oneDcenter(A[1], 1), self.oneDcenter(A[2], 2)]

    def FinetoCoarse(self, q):
        Q = [0, 0, 0]
        for x in range(3):
            if q[x] >= 0:
                Q[x] = int(self.r_0[x] + math.floor(q[x] / self.FinenodesPerCoarseCell[x]))
            else:
                Q[x] = int(self.r_0[x] + math.ceil(q[x] / self.FinenodesPerCoarseCell[x]))
        return Q

    def chooseInterpolationCoordinates(self, i, j, k, I, J, K):
        xyz = [int(np.sign(i - self.oneDcenter(I, 0))), int(np.sign(j - self.oneDcenter(J, 1))),
               int(np.sign(k - self.oneDcenter(K, 2)))]
        XYZ = [I, J, K]
        for x in range(3):
            if xyz[x] == 0:
                xyz[x] = 1
            if XYZ[x] == 0:
                xyz[x] = 1
            elif XYZ[x] == self.NumNodes[x] - 1:
                xyz[x] = -1
        Q0 = [I, J, K]
        Qx = [I + xyz[0], J, K]
        Qy = [I, J + xyz[1], K]
        Qxy = [I + xyz[0], J + xyz[1], K]
        Qz = [I, J, K + xyz[2]]
        Qyz = [I, J + xyz[1], K + xyz[2]]
        Qxz = [I + xyz[0], J, K + xyz[2]]
        Qxyz = [I + xyz[0], J + xyz[1], K + xyz[2]]
        zero = [0, 0, 0]
        if self.dim == 2:
            Q = [[Q0, Qx], [Qy, Qxy]]
            fQ = [[zero, zero], [zero, zero]]
        elif self.dim == 1:
            Q = [Q0, Qx]
            fQ = [zero, zero]
        else:
            Q = [[[Q0, Qx], [Qy, Qxy]], [[Qz, Qxz], [Qyz, Qxyz]]]
            fQ = [[[zero, zero], [zero, zero]], [[zero, zero], [zero, zero]]]
        return Q, fQ

    def linearInt(self, Q, fQ, x, axis=0):
        A, B = self.center(Q[0]), self.center(Q[1])
        fA, fB = fQ[0], fQ[1]
        fx = [0, 0, 0]
        for coord in range(3):
            fx[coord] = fA[coord] + (fB[coord] - fA[coord]) * (x[axis] - A[axis]) / (B[axis] - A[axis])
        newQ = [0, 0, 0]
        for i in range(3):
            if i == axis:
                newQ[i] = self.FinetoCoarse(x)[i]
            else:
                newQ[i] = Q[0][i]
        return newQ, fx

    def bilinearInt(self, Q, fQ, x, axis0=0, axis1=1):
        newQ1, f1 = self.linearInt(Q[0], fQ[0], x, axis0)
        newQ2, f2 = self.linearInt(Q[1], fQ[1], x, axis0)
        newQ, fx = self.linearInt([newQ1, newQ2], [f1, f2], x, axis1)
        return newQ, fx

    def quadrilinearInt(self, Q, fQ, x, axis0=0, axis1=1, axis2=2):
        newQ1, f1 = self.bilinearInt(Q[0], fQ[0], x, axis0, axis1)
        newQ2, f2 = self.bilinearInt(Q[1], fQ[1], x, axis0, axis1)
        newQ, fx = self.linearInt([newQ1, newQ2], [f1, f2], x, axis2)
        return newQ, fx

    def Interpolate(self, i, j, k, I, J, K, Field, scale=1):
        Q, fQ = self.chooseInterpolationCoordinates(i, j, k, I, J, K)
        if self.dim == 3:
            for l in range(2):
                for m in range(2):
                    for n in range(2):
                        Q_i = Q[l][m][n]
                        fQ[l][m][n] = [(Field.get(*Q_i)[0]) * scale, (Field.get(*Q_i)[1]) * scale,
                                       (Field.get(*Q_i)[2]) * scale]
            newQ, fx = self.quadrilinearInt(Q, fQ, [i, j, k])
        elif self.dim == 2:
            for l in range(2):
                for m in range(2):
                    Q_i = Q[l][m]
                    fQ[l][m] = [(Field.get(*Q_i)[0]) * scale, (Field.get(*Q_i)[1]) * scale,
                                (Field.get(*Q_i)[2]) * scale]
            newQ, fx = self.bilinearInt(Q, fQ, [i, j, k])
        elif self.dim == 1:
            for l in range(2):
                Q_i = Q[l]
                fQ[l] = [(Field.get(*Q_i)[0]) * scale, (Field.get(*Q_i)[1]) * scale, (Field.get(*Q_i)[2]) * scale]
            newQ, fx = self.linearInt(Q, fQ, [i, j, k])
        return fx

    def mixedInterpolate(self, i, j, k, I, J, K, shift, fineshift, coarseField, fineField, scale=1):
        Q, fQ = self.chooseInterpolationCoordinates(i, j, k, I, J, K)
        if self.dim == 3:
            for l in range(2):
                for m in range(2):
                    for n in range(2):
                        Q_i = Q[l][m][n]
                        if (self.r_0[0] <= Q_i[0] + shift[0] <= self.r_1[0]) and (
                                self.r_0[1] <= Q_i[1] + shift[1] <= self.r_1[1]):
                            int_i, int_j, int_k = self.center(Q_i)
                            fQ[l][m][n] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k)
                        else:
                            fQ[l][m][n] = [(coarseField.get(*Q_i)[0]) * scale, (coarseField.get(*Q_i)[1]) * scale,
                                           (coarseField.get(*Q_i)[2]) * scale]
            newQ, fx = self.quadrilinearInt(Q, fQ, [i, j, k])
        elif self.dim == 2:
            for l in range(2):
                for m in range(2):
                    Q_i = Q[l][m]
                    if (self.r_0[0] <= Q_i[0] + shift[0] <= self.r_1[0]) and (
                            self.r_0[1] <= Q_i[1] + shift[1] <= self.r_1[1]):
                        int_i, int_j, int_k = self.center(Q_i)
                        fQ[l][m] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k)
                    else:
                        fQ[l][m] = [(coarseField.get(*Q_i)[0]) * scale, (coarseField.get(*Q_i)[1]) * scale,
                                    (coarseField.get(*Q_i)[2]) * scale]
            newQ, fx = self.bilinearInt(Q, fQ, [i, j, k])
        elif self.dim == 1:
            for l in range(2):
                Q_i = Q[l]
                if (self.r_0[0] <= Q_i[0] + shift[0] <= self.r_1[0]) and (
                        self.r_0[1] <= Q_i[1] + shift[1] <= self.r_1[1]):
                    int_i, int_j, int_k = self.center(Q_i)
                    fQ[l] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k)
                else:
                    fQ[l] = [(coarseField.get(*Q_i)[0]) * scale, (coarseField.get(*Q_i)[1]) * scale,
                             (coarseField.get(*Q_i)[2]) * scale]
            newQ, fx = self.linearInt(Q, fQ, [i, j, k])
        return fx

    # TODO ? @jit()

    def fineStrayFast(self):
        M_dummy = VectorField(self.mesh)
        H_Stray = VectorField(self.mesh)
        M_dummy = self.solver.state.M
        M_dummy_numpy = self.solver.state.M.to_numpy()
       	if self.useMMstrayField == False:
          M_dummy_numpy[self.r_0[0]:self.r_1[0]+1, self.r_0[1]:self.r_1[1]+1,self.r_0[2]:self.r_1[2]+1,:] = 0
          M_dummy.from_numpy(M_dummy_numpy)
        frame= 1
        self.stray.calculate(M_dummy, H_Stray)
        stray_bc = H_Stray.to_numpy()   # bc means big coarse
        stray_mc = stray_bc[self.r_0[0]:self.r_1[0]+1,self.r_0[1]:self.r_1[1]+1, self.r_0[2]:self.r_1[2]+1]
        stray_mf = strafo.resize(stray_mc,(self.fineNodes[0],self.fineNodes[1], self.fineNodes[2], 3), order=1, mode="constant",  preserve_range=False) # mf means middle (plus frame for ghost cell interpolation)
        # now correct all areas and edges, beacuse interpolation fails at these points
        framex=int( self.FinenodesPerCoarseCell[0]/2)
        framey=int( self.FinenodesPerCoarseCell[1]/2)
        framez=int( self.FinenodesPerCoarseCell[2]/2)
        #print(stray_mc.shape)
        #print(stray_mf.shape)
        #print(framez)
        nx,ny,nz = self.fineNodes
        stray_mf[:,:,0:framez,:] = stray_mf[:,:,framez,:].reshape((nx,ny,1,3))
        stray_mf[:,:,-(framez+1):,:] = stray_mf[:,:,-(framez+1),:].reshape((nx,ny,1,3))
        stray_mf[:,0:framey,:,:] = stray_mf[:,framey,:,:].reshape((nx,1,nz,3))
        stray_mf[:,-(framey+1):,:,:] = stray_mf[:,-(framey+1),:,:].reshape((nx,1,nz,3))
        stray_mf[0:framex,:,:,:] = stray_mf[framex,:,:,:].reshape((1,ny,nz,3))
        stray_mf[-(framex+1):,:,:,:] = stray_mf[-(framex+1),:,:,:].reshape((1,ny,nz,3))
        stray_mf[0:framex,0:framey,0:framez,:]=stray_mc[0,0,0,:]
        stray_mf[0:framex,0:framey,-(framez+1):,:]=stray_mc[0,0,-1,:]

        stray_mf[0:framex,-(1+framey):,0:framez,:]=stray_mc[framex, -(1+framey), 0,:]

        stray_mf[0:framex,-(1+framey):,-(framez+1):,:]=stray_mc[framex, -(1+framey),-1,:]
        stray_mf[-(1+framex):   ,0:framey,0:framez,:]=stray_mc[-(1+framex), framey, 0,:]
        stray_mf[-(1+framex):   ,0:framey,-(framez+1):,:]=stray_mc[-(1+framex), framey,-1,:]
        stray_mf[-(1+framex):   ,-(1+framey):,0:framez,:]=stray_mc[-(1+framex), -(1+framey),0,:]
        stray_mf[-(1+framex):   ,-(1+framey):,-(framez+1):,:]=stray_mc[-(1+framex), -(1+framey),-1,:]

        #stray_mf[0:framex,-(1+framey):,0:framez,:]=stray_mc[framex, -(1+framey), framez,:]
        #stray_mf[0:framex,-(1+framey):,-(framez+1):,:]=stray_mc[framex, -(1+framey),-(1+framez),:]
        #stray_mf[-(1+framex):   ,0:framey,0:framez,:]=stray_mc[-(1+framex), framey, framez,:]
        #stray_mf[-(1+framex):   ,0:framey,-(framez+1):,:]=stray_mc[-(1+framex), framey,-(1+framez),:]
        #stray_mf[-(1+framex):   ,-(1+framey):,0:framez,:]=stray_mc[-(1+framex), -(1+framey), framez,:]
        #stray_mf[-(1+framex):   ,-(1+framey):,-(framez+1):,:]=stray_mc[-(1+framex), -(1+framey),-(1+framez),:]
        stray_to_solver =  stray_mf
        maxval0= np.amax(stray_bc[:,:,:,2] )
        minval0= np.amin(stray_bc[:,:,:,2] )
        maxval=np.amax((maxval0, minval0))
        minval=np.amin((maxval0, minval0))
        self.fineSolver.state.fineStrayField.from_numpy(stray_to_solver)


    def fineStray(self):
        startstray1 = time.time()
        M_dummy = VectorField(self.mesh)
        H_Stray = VectorField(self.mesh)
        H_Stray.fill((0, 0, 0))
        startstray2 = time.time()
        for I in range(self.NumNodes[0]):
            for J in range(self.NumNodes[1]):
                for K in range(self.NumNodes[2]):
                    if self.r_0[0] <= I <= self.r_1[0] and self.r_0[1] <= J <= self.r_1[1] and self.r_0[2] <= K <= \
                            self.r_1[2]:
                        M_dummy.set(I, J, K, (0, 0, 0))
                    else:
                        M_dummy.set(I, J, K, (self.solver.state.M.get(I, J, K)))
        startstray3 = time.time()
        self.stray.calculate(M_dummy, H_Stray)
        startstray4 =time.time()
        for i in range(self.fineNodes[0]):
            I = int(i / self.FinenodesPerCoarseCell[0])
            for j in range(self.fineNodes[1]):
                J = int(j / self.FinenodesPerCoarseCell[1])
                for k in range(self.fineNodes[2]):
                    K =  int(k/self.FinenodesPerCoarseCell[2])

                    #I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.Interpolate(i, j, k, I, J, K, H_Stray)
                    self.fineSolver.state.fineStrayField.set(i, j, k, (fx[0], fx[1], fx[2]))




                    #fx_Field.set(i,j,k ,(self.Interpolate(i, j, k, I, J, K, H_Stray)[0],self.Interpolate(i, j, k, I, J, K, H_Stray)[1],self.Interpolate(i, j, k, I, J, K, H_Stray)[2]))

        #for i in range(self.fineNodes[0]):
        #    for j in range(self.fineNodes[1]):
        #        for k in range(self.fineNodes[2]):
        #            self.fineSolver.state.fineStrayField.set(i, j, k, fx_Field.get(i,j,k))


    #					H = [0,0,0]
    #					for I in [self.r_0[0]-2, self.r_0[0]-1, self.r_1[0]+1, self.r_1[0]+2]:
    #						if 0 <= I < self.NumNodes[0]:
    #							for J in range(self.r_0[1], self.r_1[1]+1):
    #								for K in range(self.r_0[2], self.r_1[2]+1):
    #									#print I,J,K, 'fine_wall_1'
    #									L = self.Delta
    #									xyz = [(self.r_0[0] + i/self.FinenodesPerCoarseCell[0] - I)*L[0], (self.r_0[1] + j/self.FinenodesPerCoarseCell[1] - J)*L[1], (self.r_0[2] + k/self.FinenodesPerCoarseCell[2] - K)*L[2]]
    #									Mx, My, Mz  = self.solver.state.M.get(I,J,K)[0], self.fineSolver.state.M.get(i,j,k)[1], self.fineSolver.state.M.get(i,j,k)[2]
    #									Hx = magneto.Hhb(xyz[0], xyz[1], xyz[2], L[0], L[1], L[2], Mx, My, Mz)
    #									Hy = magneto.Hhb(xyz[1], xyz[2], xyz[0], L[1], L[2], L[0], Mx, My, Mz)
    #									Hz = magneto.Hhb(xyz[2], xyz[0], xyz[1], L[2], L[0], L[1], Mx, My, Mz)
    #									H[0] += Hx
    #									H[1] += Hy
    #									H[2] += Hz
    #					for J in [self.r_0[1]-2, self.r_0[1]-1, self.r_1[1]+1, self.r_1[1]+2]:
    #						if 0 <= J < self.NumNodes[1]:
    #							for I in range(self.r_0[0], self.r_1[0]+1):
    #								for K in range(self.r_0[2], self.r_1[2]+1):
    #									#print I,J,K, 'fine_wall_2'
    #									L = self.Delta
    #									xyz = [(self.r_0[0] + i/self.FinenodesPerCoarseCell[0] - I)*L[0], (self.r_0[1] + j/self.FinenodesPerCoarseCell[1] - J)*L[1], (self.r_0[2] + k/self.FinenodesPerCoarseCell[2] - K)*L[2]]
    #									Mx, My, Mz  = self.solver.state.M.get(I,J,K)[0], self.fineSolver.state.M.get(i,j,k)[1], self.fineSolver.state.M.get(i,j,k)[2]
    #									Hx = magneto.Hhb(xyz[0], xyz[1], xyz[2], L[0], L[1], L[2], Mx, My, Mz)
    #									Hy = magneto.Hhb(xyz[1], xyz[2], xyz[0], L[1], L[2], L[0], Mx, My, Mz)
    #									Hz = magneto.Hhb(xyz[2], xyz[0], xyz[1], L[2], L[0], L[1], Mx, My, Mz)
    #									H[0] += Hx
    #									H[1] += Hy
    #									H[2] += Hz
    #					self.fineSolver.state.fineStrayField.set(i,j,k, (fx[0] + H[0], fx[1] + H[1], fx[2] + H[2]))
    #
    # TODO ?@jit()

    def fineExch(self):
        left = 1
        right = 1
        top = 1
        bottom = 1
        below = 1
        above = 1
        if self.r_0[0] == 0:    left = 0
        if self.r_1[0] == self.NumNodes[0] - 1:    right = 0
        if self.r_0[1] == 0:    bottom = 0
        if self.r_1[1] == self.NumNodes[1] - 1:    top = 0
        if self.r_0[2] == 0:    below = 0
        if self.r_1[2] == self.NumNodes[2] - 1:    above = 0
        self.dummymesh = RectangularMesh(
            (self.fineNodes[0] + left + right, self.fineNodes[1] + top + bottom, self.fineNodes[2] + below + above),
            (self.l[0], self.l[1], self.l[2]))
        self.M_dummy = VectorField(self.dummymesh)
        self.M_dummy.fill((0, 0, 0))
        self.OtI_H_Exch = VectorField(self.dummymesh)
        self.OtI_H_Exch.fill((0, 0, 0))
        if self.DMI == True:
            self.OtI_H_DMI = VectorField(self.dummymesh)
            self.OtI_H_DMI.fill((0, 0, 0))
            dummyDx = VectorField(self.dummymesh)
            dummyDy = VectorField(self.dummymesh)
            dummyDz = VectorField(self.dummymesh)
            dummyDx.fill(self.fineSolver.state.Dx.get(0, 0, 0))
            dummyDy.fill(self.fineSolver.state.Dy.get(0, 0, 0))
            dummyDz.fill(self.fineSolver.state.Dz.get(0, 0, 0))
        self.mu_dummy = Field(self.dummymesh)
        self.mu_dummy.fill(self.mu)
        self.J_dummy = Field(self.dummymesh)
        self.J_dummy.fill(self.J)
        if left == 1:
            i = -1
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    # print 'M.x, ' , str(self.fineSolver.state.M.get(i,j,k)[0])+', ',str(self.fineSolver.state.M.get(i,j,k)[1]),', ',str(self.fineSolver.state.M.get(i,j,k)[2])
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
                    # else:
                    #	print 'VALUE ZERO left 1, ',str(i),', ',str(j), ', ',str(k)
                    #	time.sleep(1)

        if right == 1:
            i = self.fineNodes[0]
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
                    # else:
                    #	print "VALUE ZERO right 1"
                    #	time.sleep(1)
        if bottom == 1:
            j = -1
            for i in range(self.fineNodes[0]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    l, m, n = i, j - bottom, k
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
                    # else:
                    #	print "VALUE ZERO bottom 1"
                    #	time.sleep(1)
        if top == 1:
            j = self.fineNodes[1]
            for i in range(self.fineNodes[0]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
                    # else:
                    #	print "VALUE ZERO top 1"
                    #	time.sleep(1)
        if below == 1:
            k = -1
            for i in range(self.fineNodes[0]):
                for j in range(self.fineNodes[1]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
        #			else:
        #				print "VALUE ZERO below 1"
        #				time.sleep(1)
        if above == 1:
            k = self.fineNodes[2]
            for i in range(self.fineNodes[0]):
                for j in range(self.fineNodes[1]):
                    #			if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.mixedInterpolate(i, j, k, I, J, K, [0, 0, 0], [0, 0, 0], self.solver.state.M,
                                               self.fineSolver.state.M, scale=self.mu / self.Ms)
                    self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
        #					else:
        #						print "VALUE ZERO above 1"
        #						time.sleep(1)
        magneto.fs_exchange(self.mu_dummy, self.J_dummy, self.M_dummy, self.OtI_H_Exch)
        if self.DMI == True:
            magneto.fs_dmi(self.mu_dummy, dummyDx, dummyDy, dummyDz, self.M_dummy, self.OtI_H_DMI)

        if left == 1:
            i = 0
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                        # print 'left: Exch=', self.fineSolver.state.fineExchField.get(i, j, k), 'DMI =', self.fineSolver.state.fineDMIField.get(i, j, k)
        #			else:
        #				print "VALUE ZERO left 2"
        #				time.sleep(1)
        if right == 1:
            i = self.fineNodes[0] - 1
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    #			if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                        # print 'right: Exch=', self.fineSolver.state.fineExchField.get(i, j, k), 'DMI =', self.fineSolver.state.fineDMIField.get(i, j, k)
                    # else:
                    #	print "VALUE ZERO right 2"
                    #	time.sleep(1)

        if bottom == 1:
            j = 0
            for i in range(self.fineNodes[0]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                    # else:
                    #	print "VALUE ZERO bottom 2"
                    #	time.sleep(1)
        if top == 1:
            j = self.fineNodes[1] - 1
            for i in range(self.fineNodes[0]):
                for k in range(self.fineNodes[2]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                    # else:
                    #	print "VALUE ZERO top 2"
                    #	time.sleep(1)
        if below == 1:
            k = 0
            for i in range(self.fineNodes[0]):
                for j in range(self.fineNodes[1]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                    # else:
                    #	print "VALUE ZERO below 2"
                    #	time.sleep(1)
        if above == 1:
            k = self.fineNodes[2] - 1
            for i in range(self.fineNodes[0]):
                for j in range(self.fineNodes[1]):
                    # if(self.fineSolver.state.M.get(i,j,k)[0]!=0 or self.fineSolver.state.M.get(i,j,k)[1]!=0 or self.fineSolver.state.M.get(i,j,k)[2]!=0):
                    self.fineSolver.state.fineExchField.set(i, j, k,
                                                            self.OtI_H_Exch.get(i + left, j + bottom, k + below))
                    #if False:  # self.DMI == True:
                    if self.DMI == True:  # self.DMI == True:
                        self.fineSolver.state.fineDMIField.set(i, j, k,
                                                               self.OtI_H_DMI.get(i + left, j + bottom, k + below))
                    else:
                    	print("VALUE ZERO above 2")
                    #	time.sleep(1)
    def cellAverageFast(self,field, scale=1):
        #cellVolume=self.FinenodesPerCoarseCell[0]*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[2]
        M = VectorField(self.mesh)
        #print(scale)
        fieldnumpy =field.to_numpy()
        #print("fieldnumpy shape "+str(fieldnumpy.shape))
        avgfield = scale*strafo.downscale_local_mean(fieldnumpy,(self.FinenodesPerCoarseCell[0], self.FinenodesPerCoarseCell[1], self.FinenodesPerCoarseCell[2],1))
        #print("avgfield shape "+str(avgfield.shape))
        #print("avg field[2,2,0] "+str(avgfield[2,2,0]))
        return avgfield

    def cellAverage(self, field, scale=1):
        cellVolume = self.FinenodesPerCoarseCell[0] * self.FinenodesPerCoarseCell[1] * self.FinenodesPerCoarseCell[2]
        M = VectorField(self.mesh)
        M.fill((0, 0, 0))
        for I in range(self.r_0[0], self.r_1[0] + 1):
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    mag = [0, 0, 0]
                    for l in range(self.FinenodesPerCoarseCell[0]):
                        i = (I - self.r_0[0]) * self.FinenodesPerCoarseCell[0] + l
                        for m in range(self.FinenodesPerCoarseCell[1]):
                            j = (J - self.r_0[1]) * self.FinenodesPerCoarseCell[1] + m
                            for n in range(self.FinenodesPerCoarseCell[2]):
                                k = (K - self.r_0[2]) * self.FinenodesPerCoarseCell[2] + n
                                for x in range(3):
                                    mag[x] += (field.get(i, j, k)[x]) * (scale / cellVolume)
                    M.set(I, J, K, (mag[0], mag[1], mag[2]))
        return M

    # TODO @jit ?
    def coarseStray(self):
        M = self.cellAverage(self.fineSolver.state.M, self.Ms / self.mu)
        self.stray.calculate(M, self.solver.state.coarseStrayField)

    def coarseExch(self):
        dx = 1  # (((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**-2) -1
        dy = 1  # (((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**-2) -1
        dz = 1  # (((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**-2) -1
        scale = self.mu / self.body.material.Ms
        H_tot = VectorField(self.mesh)
        H_tot.fill((0, 0, 0))
        H_totDMI = VectorField(self.mesh)
        H_totDMI.fill((0, 0, 0))

        if self.r_0[1] != 0:
            for I in range(self.r_0[0], self.r_1[0] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    wallmesh = RectangularMesh((self.zoomsize[0] + 1, 2, self.zoomsize[2] + 1),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    if self.DMI == True:
                        H_DMI = VectorField(wallmesh)
                        H_DMI.fill((0, 0, 0))
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for i in range(self.FinenodesPerCoarseCell[0]):
                        for k in range(self.FinenodesPerCoarseCell[2]):
                            for x in range(3):
                                bottomside[x] += \
                                self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0] * (I - self.r_0[0]) + i, 0,
                                                            self.FinenodesPerCoarseCell[2] * (K - self.r_0[2]) + k)[x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[0] * self.FinenodesPerCoarseCell[
                                2])  # *((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**2)
                    oldbottom = self.solver.state.M.get(I, self.r_0[1], K)
                    M_dummy.set(I - self.r_0[0], 0, K - self.r_0[2], (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for I in range(self.r_0[0], self.r_1[0] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    H_tot.set(I, self.r_0[1] - 1, K, H_Exch.get(I - self.r_0[0], 1, K - self.r_0[2]))
                    if self.DMI == True:
                        H_totDMI.set(I, self.r_0[1] - 1, K, H_DMI.get(I - self.r_0[0], 1, K - self.r_0[2]))

        if self.r_1[1] != (self.NumNodes[1] - 1):
            for I in range(self.r_0[0], self.r_1[0] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    wallmesh = RectangularMesh((self.zoomsize[0] + 1, 2, self.zoomsize[2] + 1),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    if self.DMI == True:
                        H_DMI = VectorField(wallmesh)
                        H_DMI.fill((0, 0, 0))
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for i in range(self.FinenodesPerCoarseCell[0]):
                        for k in range(self.FinenodesPerCoarseCell[2]):
                            for x in range(3):
                                bottomside[x] += \
                                self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0] * (I - self.r_0[0]) + i,
                                                            self.fineNodes[1] - 1,
                                                            self.FinenodesPerCoarseCell[2] * (K - self.r_0[2]) + k)[x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[0] * self.FinenodesPerCoarseCell[
                                2])  # *((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**2)
                    oldbottom = self.solver.state.M.get(I, self.r_1[1], K)
                    M_dummy.set(I - self.r_0[0], 0, K - self.r_0[2], (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for I in range(self.r_0[0], self.r_1[0] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    H_tot.set(I, self.r_1[1] + 1, K, H_Exch.get(I - self.r_0[0], 1, K - self.r_0[2]))
                    if self.DMI == True:
                        H_totDMI.set(I, self.r_1[1] + 1, K, H_DMI.get(I - self.r_0[0], 1, K - self.r_0[2]))

        if self.r_0[0] != 0:
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    wallmesh = RectangularMesh((2, self.zoomsize[1] + 1, self.zoomsize[2] + 1),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    if self.DMI == True:
                        H_DMI = VectorField(wallmesh)
                        H_DMI.fill((0, 0, 0))
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for j in range(self.FinenodesPerCoarseCell[1]):
                        for k in range(self.FinenodesPerCoarseCell[2]):
                            for x in range(3):
                                bottomside[x] += \
                                self.fineSolver.state.M.get(0, self.FinenodesPerCoarseCell[1] * (J - self.r_0[1]) + j,
                                                            self.FinenodesPerCoarseCell[2] * (K - self.r_0[2]) + k)[x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[1] * self.FinenodesPerCoarseCell[
                                2])  # *((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**2)
                    oldbottom = self.solver.state.M.get(self.r_0[0], J, K)
                    M_dummy.set(0, J - self.r_0[1], K - self.r_0[2], (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, 1, 1, 1)
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    H_tot.set(self.r_0[0] - 1, J, K, H_Exch.get(1, J - self.r_0[1], K - self.r_0[2]))
                    if self.DMI == True:
                        H_totDMI.set(self.r_0[0] - 1, J, K, H_DMI.get(1, J - self.r_0[1], K - self.r_0[2]))

        if self.r_1[0] != (self.NumNodes[0] - 1):
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    wallmesh = RectangularMesh((2, self.zoomsize[1] + 1, self.zoomsize[2] + 1),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    H_DMI = VectorField(wallmesh)
                    H_DMI.fill((0, 0, 0))
                    if self.DMI == True:
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for j in range(self.FinenodesPerCoarseCell[1]):
                        for k in range(self.FinenodesPerCoarseCell[2]):
                            for x in range(3):
                                bottomside[x] += self.fineSolver.state.M.get(self.fineNodes[0] - 1,
                                                                             self.FinenodesPerCoarseCell[1] * (
                                                                                         J - self.r_0[1]) + j,
                                                                             self.FinenodesPerCoarseCell[2] * (
                                                                                         K - self.r_0[2]) + k)[x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[1] * self.FinenodesPerCoarseCell[
                                2])  # *((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**2)  TEST
                    oldbottom = self.solver.state.M.get(self.r_1[0], J, K)
                    M_dummy.set(0, J - self.r_0[1], K - self.r_0[2], (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, 1, 1, 1)  # TEST
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for K in range(self.r_0[2], self.r_1[2] + 1):
                    H_tot.set(self.r_1[0] + 1, J, K, H_Exch.get(1, J - self.r_0[1], K - self.r_0[2]))
                    if self.DMI == True:
                        H_totDMI.set(self.r_1[0] + 1, J, K, H_DMI.get(1, J - self.r_0[1], K - self.r_0[2]))

        if self.r_0[2] != 0:
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for I in range(self.r_0[2], self.r_1[2] + 1):
                    wallmesh = RectangularMesh((self.zoomsize[0] + 1, self.zoomsize[1] + 1, 2),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    if self.DMI == True:
                        H_DMI = VectorField(wallmesh)
                        H_DMI.fill((0, 0, 0))
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for i in range(self.FinenodesPerCoarseCell[0]):
                        for j in range(self.FinenodesPerCoarseCell[1]):
                            for x in range(3):
                                bottomside[x] += \
                                self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0] * (I - self.r_0[0]) + i,
                                                            self.FinenodesPerCoarseCell[1] * (J - self.r_0[1]) + j, 0)[
                                    x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[1] * self.FinenodesPerCoarseCell[
                                0])  # *((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**2)
                    oldbottom = self.solver.state.M.get(I, J, self.r_0[2])
                    M_dummy.set(I - self.r_0[0], J - self.r_0[1], 0, (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for I in range(self.r_0[0], self.r_1[0] + 1):
                    H_tot.set(I, J, self.r_0[2] - 1, H_Exch.get(I - self.r_0[0], J - self.r_0[1], 1))
                    if self.DMI == True:
                        H_totDMI.set(I, J, self.r_0[2] - 1, H_DMI.get(I - self.r_0[0], J - self.r_0[1], 1))

        if self.r_1[2] != (self.NumNodes[2] - 1):
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for I in range(self.r_0[0], self.r_1[0] + 1):
                    wallmesh = RectangularMesh((self.zoomsize[0] + 1, self.zoomsize[1] + 1, 2),
                                               (self.Delta[0], self.Delta[1], self.Delta[2]))
                    M_dummy = VectorField(wallmesh)
                    M_dummy.fill((0, 0, 0))
                    Ms_dummy = Field(wallmesh)
                    Ms_dummy.fill(self.body.material.Ms)
                    A_dummy = Field(wallmesh)
                    A_dummy.fill(self.body.material.A)
                    H_Exch = VectorField(wallmesh)
                    H_Exch.fill((0, 0, 0))
                    if self.DMI == True:
                        H_DMI = VectorField(wallmesh)
                        H_DMI.fill((0, 0, 0))
                        wallDx = VectorField(wallmesh)
                        wallDy = VectorField(wallmesh)
                        wallDz = VectorField(wallmesh)
                        wallDx.fill(self.solver.state.Dx.get(0, 0, 0))
                        wallDy.fill(self.solver.state.Dy.get(0, 0, 0))
                        wallDz.fill(self.solver.state.Dz.get(0, 0, 0))
                    bottomside = [0, 0, 0]
                    bottomrescaled = [0, 0, 0]
                    for i in range(self.FinenodesPerCoarseCell[0]):
                        for k in range(self.FinenodesPerCoarseCell[2]):
                            for x in range(3):
                                bottomside[x] += \
                                self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0] * (I - self.r_0[0]) + i,
                                                            self.FinenodesPerCoarseCell[1] * (J - self.r_0[1]) + j,
                                                            self.fineNodes[2] - 1)[x]
                    for x in range(3):
                        bottomrescaled[x] = bottomside[x] / (
                                    scale * self.FinenodesPerCoarseCell[1] * self.FinenodesPerCoarseCell[
                                0])  # *((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**2)
                    oldbottom = self.solver.state.M.get(I, J, self.r_1[2])
                    M_dummy.set(I - self.r_0[0], J - self.r_0[1], 0, (
                    bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1],
                    bottomrescaled[2] - oldbottom[2]))
            magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
            if self.DMI == True:
                magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
            for J in range(self.r_0[1], self.r_1[1] + 1):
                for I in range(self.r_0[0], self.r_1[0] + 1):
                    H_tot.set(I, J, self.r_1[2] + 1, H_Exch.get(I - self.r_0[0], J - self.r_0[1], 1))
                    if self.DMI == True:
                        H_totDMI.set(I, J, self.r_1[2] + 1, H_DMI.get(I - self.r_0[0], J - self.r_0[1], 1))

        self.solver.state.coarseExchField = H_tot

        if self.DMI == True:
            self.solver.state.coarseDMIField = H_totDMI

    def fillFineScale(self):
        for i in range(self.fineNodes[0]):
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    fx = self.Interpolate(i, j, k, I, J, K, self.solver.state.M, self.mu / self.Ms)
                    self.fineSolver.state.M.set(i, j, k, (fx[0], fx[1], fx[2]))
        self.fineMagIsGiven = True

    def fillFineScaleFast(self):
        frame = 1
        solvernumpy = self.solver.state.M.to_numpy()
        fineSolverPlusFrame = solvernumpy[self.r_0[0]-frame:self.r_1[0]+frame+1, self.r_0[1]-frame:self.r_1[1]+frame+1, self.r_0[2]-2:self.r_1[2]+frame+1]
        fineSolverPlusFrameInterpolated = strafo.resize(fineSolverPlusFrame,(self.fineNodes[0]+(2*frame)*self.FinenodesPerCoarseCell[0],self.fineNodes[1]+(2*frame)*self.FinenodesPerCoarseCell[1], self.FinenodesPerCoarseCell[2] ,3), order = 1, mode = "symmetric", preserve_range= True)
        self.fineSolver.state.M.from_numpy(self.mu/self.Ms*fineSolverPlusFrameInterpolated[self.FinenodesPerCoarseCell[0]*frame: - frame*self.FinenodesPerCoarseCell[0],self.FinenodesPerCoarseCell[1]*frame: - frame*self.FinenodesPerCoarseCell[1],self.FinenodesPerCoarseCell[2]])

        self.fineMagIsGiven = True

    def Multiscaleloop(self):
        self.h = self.solver.state.h
        self.t = self.solver.state.t
        self.fineSolver.state.fineExchField.fill((0, 0, 0))
        if self.DMI == True:
            self.fineSolver.state.fineDMIField.fill((0, 0, 0))
        self.startfinestray = time.time()
        self.fineStrayFast()
        self.stopfinestray = time.time()
        fexstart=time.time()
        self.fineExch()
        fexstop=time.time()

        if not self.t == 0:
            #if self.multiscaleLog:
                #print("_____________________________  coarse/fine step: "+str(self.solver.state.step)+","+str(self.fineSolver.state.step))
            self.startsolve=time.time()
            #if self.fineSolver.state.step != 1:
            #    self.fineSolver.solve(condition.Time(self.t + (self.h) / 2))
            #else:
            #print("stateM(avg): "+str(self.fineSolver.state.M.get(0,0,0)))

            #print("stateM(avg): "+str(self.fineSolver.state.M.get(310,310,0)))
            #print("stateM(avg): "+str(self.fineSolver.state.M.get(10,10,0)))
            self.fineSolver.solve(condition.Time(self.t+(self.h)/2))
            #print("solver stray 200 200", self.solver.state.H_stray.get(200,200,0))
            #print("solver stray 200 201", self.solver.state.H_stray.get(200,201,0))

            #print("finesolver stray 202 202", self.fineSolver.state.fineStrayField.get(202,202,0))

            #print("finesolver stray 202 203", self.fineSolver.state.fineStrayField.get(202,203,0))

            #self.fineSolver.solve(condition.Time(self.t + (self.h) / 2))

            #print("state after: "+str(self.fineSolver.state.M.get(0,0,0)))

            #print("stateM(avg): "+str(self.fineSolver.state.M.get(310,310,0)))
            #print("stateM(avg): "+str(self.fineSolver.state.M.get(10,10,0)))
            self.stopsolve = time.time()
            # self.solver.state.fine_deg_per_ns = self.fineSolver.state.deg_per_ns
            setattr(self.solver.state, "fine_deg_per_ns", self.fineSolver.state.deg_per_ns)
            # print(self.solver.state.__dict__)
        else:
            self.f = open(self.folder + '/pos.txt', 'w')
            self.f.write(str([self.t, self.r_0, self.r_1]) + '\n')
        #Mavg = self.cellAverage(self.fineSolver.state.M, self.Ms / self.mu)

        #for I in range(self.NumNodes[0]):
        #    for J in range(self.NumNodes[1]):
        #        for K in range(self.NumNodes[2]):
        #            if not (self.r_0[0] <= I <= self.r_1[0] and self.r_0[1] <= J <= self.r_1[1] and self.r_0[2] <= K <=
        #                    self.r_1[2]):
        #                Mavg.set(I, J, K, (self.solver.state.M.get(I, J, K)))
        #self.solver.state.M = Mavg
        #self.solver.state.coarseStrayField.fill((0, 0, 0))
        #self.solver.state.coarseExchField.fill((0, 0, 0))
        #if self.DMI == True:
        #    self.solver.state.coarseDMIField.fill((0, 0, 0))
        #self.coarseExch()
        #self.done()
        #print(self.solver.state.M.get(self.r_0[0]+2,self.r_0[1]+2,self.r_0[2]+0))

        Mavg = self.cellAverageFast(self.fineSolver.state.M, self.Ms/self.mu)
        #print(np.shape(Mavg))
        #time1 = time.time()
        solverNum = self.solver.state.M.to_numpy()
        #print(np.shape(solverNum))

        solverNum[self.r_0[0]:self.r_1[0]+1,self.r_0[1] :self.r_1[1]+1,self.r_0[2] :self.r_1[2]+1] = Mavg
        self.solver.state.M.from_numpy(solverNum)
        #time2=time.time()
        #print("r_0")
        #print(self.solver.state.M.get(self.r_0[0],self.r_0[1],self.r_0[2]) )
        #for i in range(self.r_1[0]-self.r_0[0]+1):
        #    for j in range(self.r_1[1]-self.r_0[1]+1):
        #        for k in range(self.r_1[2]-self.r_0[2]+1):
                    #print(self.fineNodes[1])

         #          self.solver.state.M.set(self.r_0[0]+i,self.r_0[1]+j,self.r_0[2]+k,(Mavg[i, j, k,0],Mavg[i, j, k,1],Mavg[i, j, k,2]))

        #print(self.solver.state.M.get(self.r_0[0],self.r_0[1],self.r_0[2] ))
        #time3 = time.time()
        #print("fast: "+str(time2-time1))
        #print("slow: "+str(time3-time2))
        #avg = self
        self.solver.state.coarseStrayField.fill((0,0,0))
        self.solver.state.coarseExchField.fill((0,0,0))
        if self.DMI == True:
            self.solver.state.coarseDMIField.fill((0,0,0))
        self.coarseExch()
        time6 = time.time()
        #self.wholefineExch += fexstop -fexstart
        #print("whole coarseExch "+str(self.wholecoarseExch))
        #print("whole fineExch "+str(self.wholefineExch))
        self.done()

    def shrinkright(self, shrinkcoarserangex):
        # shrink T/B stephandler
        # copying to tmp values
        print("shrinkcoarserangex: "+str(shrinkcoarserangex))
        print("here we are, shrinkright " +self.TrackingCond)
        shrinkfinerangex = shrinkcoarserangex * self.FinenodesPerCoarseCell[0]
        tmpr1 = [self.r_1[0] - abs(shrinkcoarserangex), self.r_1[1], self.r_1[2]]
        tmpr0 = self.r_0
        tmpfineNodes = (self.fineNodes[0] - abs(shrinkfinerangex), self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (self.l[0], self.l[1], self.l[2]))
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))



        #self.fineSolver.system.modules.remove(LandauLifshitzGilbert)
        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        #del tmpmodulelist[0]
        #not with pop
        #tmpmodulelist = tmpmodulelist[1:]
        #tmpmodulelist.remove(LandauLifshitzGilbert)
                    ###remove(0# )
                    ###finesolver.module.remove(LLG) -> not in list
        #print("module list leng after pop:"+str(len(tmpmodulelist)))
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        print('mancopy: '+str(mancopy))
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        tmpfineSolver.state.t = self.fineSolver.state.t
        tmpfineSolver.state.h = self.fineSolver.state.h
        tmpfineSolver.state.step = self.fineSolver.state.step
        print("simple external field value:"+str(self.fineSolver.state.fineExternalField.get(0,0,0)))

        tmpfineSolver.state.fineExternalField.fill((self.fineSolver.state.fineExternalField.get(0,0,0)))
        print("manually copied"+str(self.fineSolver.state.fineExternalField.uniform_value))


        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfineSolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i, j, k))        # TODO tmpbuffer einfuege


        self.reInit( tmpr0, tmpr1, tmpfineSolver, mancopy)

    def growright(self, growcoarserangex):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, growight")
        growfinerangex = growcoarserangex * self.FinenodesPerCoarseCell[0]
        tmpr1 = [self.r_1[0] + growcoarserangex, self.r_1[1], self.r_1[2]]
        tmpr0 = self.r_0
        tmpFineNodes = (self.fineNodes[0] + growfinerangex, self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpFineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody=self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!

        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)
        tmpfineSolver.state.fineExternalField.fill(self.fineSolver.state.fineExternalField.uniform_value)
        for i in range(tmpFineNodes[0]):
            for j in range(tmpFineNodes[1]):
                for k in range(tmpFineNodes[2]):
                    tmpfineSolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i, j, k))

        # TODO scale!
        for i in range(self.fineNodes[0], self.fineNodes[0] + growfinerangex):
            for j in range(tmpFineNodes[1]):
                for k in range(tmpFineNodes[2]):
                    I, J, K = self.FinetoCoarse(i, j, k)
                    fx = self.Interpolate(i, j, k, I, J, K, self.solver.state.M, scale=1)
                    tmpfineSolver.state.M.set(i, j, k, fx[0], fx[1], fx[2])
        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver)

    def shrinktop(self, shrinkcoarserangey):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, shrinktop "+self.TrackingCond)
        shrinkfinerangey = abs(shrinkcoarserangey) * self.FinenodesPerCoarseCell[0]
        tmpr1 = [self.r_1[0], self.r_1[1] - abs(shrinkcoarserangey), self.r_1[2]]
        tmpr0 = self.r_0
        tmpfineNodes = (self.fineNodes[0], self.fineNodes[1] - abs(shrinkfinerangey), self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody = self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!
        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfineSolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i, j, k))

        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

    def growtop(self, growcoarserangey):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, growtop " +self.TrackingCond)
        growfinerangey = growcoarserangey * self.FinenodesPerCoarseCell[0]
        tmpr1 = [self.r_1[0], self.r_1[1] + growcoarserangey, self.r_1[2]]
        tmpr0 = self.r_0
        tmpfineNodes = (self.fineNodes[0], self.fineNodes[1] + growfinerangey, self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody = self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!

        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfinesolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i, j, k))
        # TODO scale=1?!
        for i in range(tmpfineNodes[0]):
            for j in range(self.fineNodes[1], self.fineNodes[1] + growfinerangey):
                for k in range(tmpfineNodes[2]):
                    I, J, K = self.FinetoCoarse((i, j, k))
                    fx = self.Interpolate(i, j, k, I, J, K, self.solver.state.M, scale=1)
                    tmpfineSolver.state.M.set(i, j, k, fx[0], fx[1], fx[2])

        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

        # TODO mixedInterpolate Andrea fragen!!

    def growleft(self, growcoarserangex):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, growleft "+self.TrackingCond)
        growfinerangex = abs(growcoarserangex) * self.FinenodesPerCoarseCell[0]
        tmpr1 = self.r_1
        tmpr0 = [self.r_0[0] - abs(growcoarserangex), self.r_0[1], self.r_0[2]]
        tmpfineNodes = (self.fineNodes[0] + abs(growfinerangex), self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody=self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!
        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfineSolver.state.M.set(i + growfinerangex, j, k, self.fineSolver.state.M.get(i, j, k))

        # TODO scale!
        for i in range(growfinerangex):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    I, J, K = self.FinetoCoarse((i, j, k))
                    fx = self.Interpolate(i, j, k, I, J, K, self.solver.state.M, scale=1)
                    tmpfineSolver.state.M.set((i, j, k), fx[0], fx[1], fx[2])
        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

    def shrinkleft(self, shrinkcoarserangex):

        print("here we are, shrinkleft "+self.TrackingCond)
        shrinkfinerangex = abs(shrinkcoarserangex) * self.FinenodesPerCoarseCell[0]
        tmpr1 = self.r_1
        tmpr0 = [self.r_0[0] + abs(shrinkcoarserangex), self.r_0[1], self.r_0[2]]
        tmpfineNodes = (self.fineNodes[0] - shrinkfinerangex, self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody=self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))

        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        print("mancopy: "+str(mancopy) )

        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)
        tmpfineSolver.state.t = self.fineSolver.state.t
        tmpfineSolver.state.h = self.fineSolver.state.h
        tmpfineSolver.state.step = self.fineSolver.state.step
        print("manually copied"+str(self.fineSolver.state.fineExternalField.uniform_value))

        #tmpfineSolver.evolver.__runge_kutta_next_h = self.fineSolver.evolver.__runge_kutta_next_h
        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfineSolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i + shrinkfinerangex, j, k))

        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

    def shrinkbottom(self, shrinkcoarserangey):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, shrinkbottom "+self.TrackingCond)
        shrinkfinerangey = abs(shrinkcoarserangey) * self.FinenodesPerCoarseCell[0]
        tmpr1 = self.r_1
        tmpr0 = [self.r_0[0], self.r_0[1] + abs(shrinkcoarserangey), self.r_0[2]]
        tmpfineNodes = (self.fineNodes[0] - abs(shrinkfinerangey), self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody=self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!

        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpfineSolver.state.M.set(i, j, k, self.fineSolver.state.M.get(i , j+ shrinkfinerangey, k))

        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

    def growbottom(self, growcoarserangey):
        # shrink T/B stephandler
        # copying to tmp values
        print("here we are, growbottom "+self.TrackingCond)
        growfinerangey = abs(growcoarserangey) * self.FinenodesPerCoarseCell[0]
        tmpr1 = self.r_1
        tmpr0 = [self.r_0[0], self.r_0[1] - abs(growcoarserangey), self.r_0[2]]
        tmpfineNodes = (self.fineNodes[0] + abs(growfinerangey), self.fineNodes[1], self.fineNodes[2])
        tmpmesh = RectangularMesh(tmpfineNodes, (
            self.l[0], self.l[1], self.l[2]))
        # tmpbody=self.fineSolver.body
        tmpbodyID=self.fineWorld.bodies[0]._Body__id
        TMPDISK =  self.fineWorld.bodies[0]._Body__material
        tmpworld = World(tmpmesh, Body(tmpbodyID, TMPDISK))
        # TODO lambda functions veraendern, da mesh veraendert ist!!!
        ## TODO geht auch nicht mit solver.mesh!! weil die  nicht ueber self. ansprechbar sind
        # TODO modulelist uebernehmen!

        tmpmodulelist = self.fineSolver.system.modules.__getslice__(1,len(self.fineSolver.system.modules)+1)
        mancopy = self.fineSolver.state.fineExternalField.uniform_value
        #defining mancopy must be done before create tmpsolver (don know why)
        tmpfineSolver = create_solver(tmpworld, tmpmodulelist, True, log=True)

        for i in range(tmpfineNodes[0]):
            for j in range(tmpfineNodes[1]):
                for k in range(tmpfineNodes[2]):
                    tmpsolver.state.M.set(i, j+ shrinkfinerangey, k, self.fineSolver.state.M.get(i , j, k))

        for i in range(tmpfineNodes[0]):
            for j in range(shrinkfinerangey):
                for k in range(tmpfineNodes[2]):
                    I, J, K = self.FinetoCoarse((i, j, k))
                    fx = self.Interpolate(i, j, k, I, J, K, self.f)
                    tmpsolver.state.M.set(i, j, k, (fx[0], fx[1], fx[2]))

        # TODO tmpbuffer einfuegen
        self.reInit(tmpr0, tmpr1, tmpfineSolver, mancopy)

    def shift(self):
        self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
        self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
        self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0],
                            self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]
        self.newCoarseshift = [0, 0, 0]
        for x in range(2):
            if self.r_0[x] + self.Coarseshift[x] < 0:
                self.newCoarseshift[x] = - self.r_0[x]
            elif self.r_1[x] + self.Coarseshift[x] > self.NumNodes[x] - 1:
                self.newCoarseshift[x] = self.NumNodes[x] - 1 - self.r_1[x]
            else:
                self.newCoarseshift[x] = self.Coarseshift[x]
            self.r_0[x] += self.newCoarseshift[x]
            self.r_1[x] += self.newCoarseshift[x]
            self.Fineshift[x] = self.newCoarseshift[x] * self.FinenodesPerCoarseCell[x]
        new_M = VectorField(self.fineMesh)
        for i in range(self.fineNodes[0]):
            new_i = i + self.Fineshift[0]
            for j in range(self.fineNodes[1]):
                new_j = j + self.Fineshift[1]
                for k in range(self.fineNodes[2]):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    if 0 <= new_i < self.fineNodes[0] and 0 <= new_j < self.fineNodes[1]:
                        new_M.set(i, j, k, self.fineSolver.state.M.get(new_i, new_j, k))
                    else:
                        fx = self.mixedInterpolate(i, j, k, I, J, K, self.newCoarseshift, self.Fineshift,
                                                   self.solver.state.M, self.fineSolver.state.M,
                                                   scale=self.mu / self.Ms)
                        new_M.set(i, j, k, (fx[0], fx[1], fx[2]))
        self.fineSolver.state.M = new_M
        self.f.write(str([self.t, self.r_0, self.r_1]) + '\n')
        self.oldcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.newcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
        self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
        self.Coarseshift = [0, 0, 0]
        self.Multiscaleloop()

    def shiftForced(self, shiftcoarserangex, shiftcoarserangey, shiftcoarserangez):
        ## TODO: When becomes newcenter = colcenter??? (should be done)
        self.Coarseshift = [shiftcoarserangex, shiftcoarserangey, shiftcoarserangez]
        print("shift Forced "+ str(self.TrackingCond))
        print(self.Coarseshift)
        print(self.oldcenterCoarse)
        self.newCoarseshift = [0, 0, 0]
        for x in range(2):
            if self.r_0[x] + self.Coarseshift[x] < 0:
                self.newCoarseshift[x] = - self.r_0[x]
            elif self.r_1[x] + self.Coarseshift[x] > self.NumNodes[x] - 1:
                self.newCoarseshift[x] = self.NumNodes[x] - 1 - self.r_1[x]
            else:
                self.newCoarseshift[x] = self.Coarseshift[x]
            self.r_0[x] += self.newCoarseshift[x]
            self.r_1[x] += self.newCoarseshift[x]
            self.Fineshift[x] = self.newCoarseshift[x] * self.FinenodesPerCoarseCell[x]
        new_M = VectorField(self.fineMesh)
        for i in range(self.fineNodes[0]):
            new_i = i + self.Fineshift[0]
            for j in range(self.fineNodes[1]):
                new_j = j + self.Fineshift[1]
                for k in range(self.fineNodes[2]):
                    I, J, K = self.FinetoCoarse([i, j, k])
                    if 0 <= new_i < self.fineNodes[0] and 0 <= new_j < self.fineNodes[1]:
                        new_M.set(i, j, k, self.fineSolver.state.M.get(new_i, new_j, k))
                    else:
                        fx = self.mixedInterpolate(i, j, k, I, J, K, self.newCoarseshift, self.Fineshift,
                                                   self.solver.state.M, self.fineSolver.state.M,
                                                   scale=self.mu / self.Ms)
                        new_M.set(i, j, k, (fx[0], fx[1], fx[2]))
        self.fineSolver.state.M = new_M
        self.f.write(str([self.t, self.r_0, self.r_1]) + '\n')
        # TODO: following neccessary??
        #print("shift forced " + str(self.TrackingCond))
        #print("r_0 : "+str(self.r_0) )
        self.Coarseshift = [0, 0, 0]
        self.oldcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.newcenter = [math.floor(self.fineNodes[0] / 2), math.floor(self.fineNodes[1] / 2),
                          math.floor(self.fineNodes[2] / 2)]
        self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
        self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
        #print("new center")
        #print(self.oldcenterCoarse)

    def trackVortex(self):
        mini = 0
        maxi = 0
        coordmini = [0, 0]
        coordmaxi = [0, 0]
        B = int(self.fineNodes[0] * 0.1)
        for i in range(int(self.oldcenter[0] - B), int(self.oldcenter[0] + B)):
            for j in range(int(self.oldcenter[1] - B), int(self.oldcenter[1] + B)):
                A = 0
                for k in range(self.fineNodes[2]):
                    A += self.fineSolver.state.M.get(i, j, k)[2] / (self.fineNodes[2] * self.body.material.mu)
                if A < mini:
                    mini = A
                    coordmini = [i, j, self.oldcenter[2]]
                elif A > maxi:
                    maxi = A
                    coordmaxi = [i, j, self.oldcenter[2]]
        if (mini < -0.9) and (maxi > 0.9):
            x1 = min(coordmini[0], coordmaxi[0])
            x2 = max(coordmini[0], coordmaxi[0])
            y1 = min(coordmini[1], coordmaxi[1])
            y2 = max(coordmini[1], coordmaxi[1])
            self.newcenter = [x1 + int((x2 - x1) / 2), y1 + int((y2 - y1) / 2), self.newcenter[2]]
        else:
            if maxi > abs(mini):
                self.newcenter = coordmaxi
            elif abs(mini) > maxi:
                self.newcenter = coordmini
        shift = np.sqrt((self.newcenter[0] - self.oldcenter[0]) ** 2 + (self.newcenter[1] - self.oldcenter[1]) ** 2)
        C = int(self.FinenodesPerCoarseCell[0])
        if not (C < shift < B):
            self.newcenter = self.oldcenter
            self.Multiscaleloop()
        else:
            self.shift()
            self.Multiscaleloop()

    def trackDW(self):
        mini = 1
        B = int(self.fineNodes[0] * 0.1)
        for i in range(int(self.oldcenter[0] - B), int(self.oldcenter[0] + B)):
            A = 0
            for j in range(self.fineNodes[1]):
                for k in range(self.fineNodes[2]):
                    A += abs(self.fineSolver.state.M.get(i, j, k)[0]) / (
                                self.fineNodes[1] * self.fineNodes[2] * self.body.material.mu)
            if A < mini:
                mini = A
                self.newcenter[0] = i
        C = int(self.fineNodes[0] * self.mintrack)
        shift = abs(self.newcenter[0] - self.oldcenter[0])
        if not (C < shift < B):
            self.newcenter[0] = self.oldcenter[0]
            self.Multiscaleloop()
        else:
            if self.t == 0:
                self.Multiscaleloop()
            else:
                self.shift()
                # self.Multiscaleloop()

    def trackSkyrmionEdgeL(self):
        # x-Direction
        # take care that x-range is smaller than y range!

        # DWin = False  # checks if the SKyrmion DW is inside the finesolver (to avoid arbitrary
        # isactive = False  # TODO if DW completely out (in passive solver) => change condition to Every 10th step!
        lowMz = 1
        #innerMz = 0
        lowMzbound = 0.04  # TODO check that, optimize value
        #innerMzbound = 0.9  # TODO check that, optimize value

        for i in range(self.fineNodes[1] - 1):
            tmpMz = abs(self.fineSolver.state.M.get(i, int(self.oldcenter[1]), 0)[2])
            lowMz = min(lowMz, tmpMz)
        lowMz /= self.mu
        # TODO check wether this makes sense (set inactive)
        #innerMz = abs(self.fineSolver.state.M.get(self.fineNodes[0] - 1, int(self.oldcenter[1]), 0)[2]) / self.mu
        #print("lowMz:" + str(lowMz))
        #print("innerMz : " + str(innerMz))

            # TODO where change condition?
            # TODO when inactive --> when active again set finescale again from coarse, set time again!!

        # find 1d-center position
        mini = 1
        Bx = int(self.fineNodes[0] * 0.3)  # to have a square area that gets evaluated
        By = int(self.fineNodes[1] * 0.04)

        for i in range(int(self.oldcenter[0]) - Bx, int(self.oldcenter[0]) + Bx, 1):
            A = 0
            for j in range(int(self.oldcenter[1]) - By, int(self.oldcenter[1]) + By, 1):
                for k in range(self.fineNodes[2]):
                    A += abs(self.fineSolver.state.M.get(i, j, k)[2]) / (
                                self.fineNodes[1] * self.fineNodes[2] * self.body.material.mu)
            if A < mini:
                mini = A
                self.newcenter[0] = i + 4
        C = int(self.fineNodes[0] * self.mintrack)
        shift = abs(self.newcenter[0] - self.oldcenter[0])
        if (lowMz < lowMzbound):

            if not (C < shift < self.fineNodes[0]):
                self.newcenter[0] = self.oldcenter[0]


            else:
                print('shift left to right from... to ...')
                print(self.oldcenter[0])
                print(i)
                print("C = fineNodes*mintrack")
                print(C)
                print('shift left to right')
                self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
                self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
                self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0],
                                    self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]
                print("Coars shif")
                print(self.Coarseshift)
                print("coarse shift corrected")
                if self.Coarseshift[0] > 1:
                    self.Coarseshift[0] = 1
                    self.newcenterCoarse[0] = self.oldcenterCoarse[0] + self.Coarseshift[0]
                if self.Coarseshift[0] < -1:
                    self.Coarseshift[0] = -1
                    self.newcenterCoarse[0] = self.oldcenterCoarse[0] + self.Coarseshift[0]
        else:

            print("DW not in finescale!")
        ## IMPORTANT: Shifting is done in the solver track Skrmion DW routing!
        #      self.shiftWithoutLoop()
        # self.Multiscaleloop()

    def trackSkyrmionEdgeR(self):

        # DWin = False  # checks if the SKyrmion DW is inside the finesolver (to avoid arbitrary shifting)
        lowMz = 1
        #innerMz = 0
        lowMzbound = 0.04  # TODO check that, optimize value
        #innerMzbound = 0.9  # TODO check that, optimize value

        for i in range(self.fineNodes[1] - 1):
            tmpMz = abs(self.fineSolver.state.M.get(i, int(self.oldcenter[1]), 0)[2])
            lowMz = min(lowMz, tmpMz)
        lowMz /= self.mu

        # TODO check with Kai wether this makes sense (set inactive)
        innerMz = abs(self.fineSolver.state.M.get(0, int(self.oldcenter[1]), 0)[2]) / self.mu

        #print("lowMz:" + str(lowMz))
        #print("innerMz : " + str(innerMz))

        # y-Direction
        # take care that x-range is smaller than y range!
        mini = 1
        Bx = int(self.fineNodes[0] * 0.3)  # to have a square area that gets evaluated
        By = int(self.fineNodes[1] * 0.04)
        # print("here we are now! x-direc of SEdge R ")
        for i in range(int(self.oldcenter[0]) - Bx, int(self.oldcenter[0]) + Bx, 1):
            A = 0
            for j in range(int(self.oldcenter[1]) - By, int(self.oldcenter[1]) + By, 1):
                for k in range(self.fineNodes[2]):
                    A += abs(self.fineSolver.state.M.get(i, j, k)[2]) / (
                            self.fineNodes[1] * self.fineNodes[2] * self.body.material.mu)
            if A < mini:
                mini = A
                self.newcenter[0] = i - 4
        C = int(self.fineNodes[0] * self.mintrack)
        shift = abs(self.newcenter[0] - self.oldcenter[0])
        if lowMz < lowMzbound:
            if not (C < shift < self.fineNodes[0]):
                self.newcenter[0] = self.oldcenter[0]
            else:
                self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
                self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
                self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0],
                                    self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]
                #print("Coars shif")
                #print(self.Coarseshift)
                print("coarse shift corrected")
                if self.Coarseshift[0] > 1:
                    self.Coarseshift[0] = 1
                    self.newcenterCoarse[0] = self.oldcenterCoarse[0] + self.Coarseshift[0]
                if self.Coarseshift[0] < -1:
                    self.Coarseshift[0] = -1
                    self.newcenterCoarse[0] = self.oldcenterCoarse[0] + self.Coarseshift[0]
        else:

            print("DW not in!")
            #       self.shiftWithoutLoop()

    def trackSkyrmionEdgeT(self):
        print("TRACK NOW!")
        # DWin = False  # checks if the SKyrmion DW is inside the finesolver (to avoid arbitrary
        # isactive = False  # TODO if DW completely out (in passive solver) => change condition to Every 10th step!
        lowMz = 1
        # innerMz = 0
        lowMzbound = 0.04  # TODO check that, optimize value

        for j in range(self.fineNodes[0] - 1):
            tmpMz = abs(self.fineSolver.state.M.get(int(self.oldcenter[0]), j, 0)[2])
            lowMz = min(lowMz, tmpMz)

        # TODO check with Kai wether this makes sense (set inactive)
        #innerMz = abs(self.fineSolver.state.M.get(int(self.oldcenter[0]), 0, 0)[2]) / self.mu

        #print("lowMz:" + str(lowMz))
        #print("innerMz : " + str(innerMz))


    # y-Direction
        mini = 1
        Bx = int(self.fineNodes[0] * 0.04)  # to have a square area which gets evaluated
        By = int(self.fineNodes[1] * 0.3)
        #print("By:" + str(By))
        for j in range(int(self.oldcenter[1]) - By, int(self.oldcenter[1]) + By, 1):
            A = 0
            for i in range(int(self.oldcenter[0]) - Bx, int(self.oldcenter[0]) + Bx, 1):
                for k in range(self.fineNodes[2]):
                    A += abs(self.fineSolver.state.M.get(i, j, k)[2]) / (
                            self.fineNodes[0] * self.fineNodes[2] * self.body.material.mu)
            if A < mini:
                mini = A
                self.newcenter[1] = j - 4
        C = int(self.fineNodes[1] * self.mintrack)
        shift = abs(self.newcenter[1] - self.oldcenter[1])
        if lowMz < lowMzbound:
            print("DW in!")
            if not (C < shift < self.fineNodes[1]):
                self.newcenter[1] = self.oldcenter[1]
                # self.Multiscaleloop()
            # else:
            # if self.t == 0:
            #    self.Multiscaleloop()
            else:
                #print('shift top down from... to ...')
                print(self.oldcenter[1])
                #print(j)
                #print("C = fineNodes*mintrack")
                #print(C)
                self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
                self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
                self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0],
                                    self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]
                #print("Coars shif")
                #print(self.Coarseshift)
                #print("coarse shift corrected")
                if self.Coarseshift[1] > 1:
                    self.Coarseshift[1] = 1
                    self.newcenterCoarse[1] = self.oldcenterCoarse[1] + self.Coarseshift[1]
                if self.Coarseshift[1] < -1:
                    self.Coarseshift[1] = -1
                    self.newcenterCoarse[1] = self.oldcenterCoarse[1] + self.Coarseshift[1]
                print("T Coarseshift: " +str(self.Coarseshift))
        else:
            print("DW not in!")

    def trackSkyrmionEdgeB(self):

        # DWin = False  # checks if the SKyrmion DW is inside the finesolver (to avoid arbitrary
        # isactive = False  # TODO if DW completely out (in passive solver) => change condition to Every 10th step!
        lowMz = 1
        innerMz = 0
        lowMzbound = 0.04  # TODO check that, optimize value
        innerMzbound = 0.9  # TODO check that, optimize value

        for j in range(self.fineNodes[0] - 1):
            tmpMz = abs(self.fineSolver.state.M.get(int(self.oldcenter[0]), j, 0)[2])
            lowMz = min(lowMz, tmpMz)
        lowMz /= self.mu

        # TODO check with Kai wether this makes sense (set inactive), scale mu / Ms
        #innerMz = abs(self.fineSolver.state.M.get(int(self.oldcenter[0]), self.fineNodes[1] - 1, 0)[2]) / self.mu

        #print("lowMz:" + str(lowMz))
        #print("innerMz : " + str(innerMz))

        # y-Direction
        mini = 1
        Bx = int(self.fineNodes[0] * 0.06)  # to have a square area which gets evaluated
        By = int(self.fineNodes[1] * 0.3)
        for j in range(int(self.oldcenter[1]) - By, int(self.oldcenter[1]) + Bx, 1):
            A = 0
            for i in range(int(self.oldcenter[0]) - Bx, int(self.oldcenter[0]) + Bx, 1):
                for k in range(self.fineNodes[2]):
                    A += abs(self.fineSolver.state.M.get(i, j, k)[2]) / (
                            self.fineNodes[0] * self.fineNodes[2] * self.body.material.mu)
            if A < mini:
                mini = A
                self.newcenter[1] = j +4
        C = int(self.fineNodes[1] * self.mintrack)

        signedshift = self.newcenter[1] - self.oldcenter[1]
        shift = abs(signedshift)
        self.newcenterCoarse

        if lowMz < lowMzbound:

            if not (C < shift < self.fineNodes[1]):
                self.newcenter[1] = self.oldcenter[1]
            #   self.Multiscaleloop()
            else:
                print('shift bottom up from ... to...')
                print(self.oldcenter[1])
                print(j)
                print("C = fineNodes*mintrack")
                print(C)
                print('shift bottom up')
                self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
                self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
                self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0],
                                    self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]


                print("coarse shift corrected")
                if self.Coarseshift[1] > 1:
                    self.Coarseshift[1] = 1
                    self.newcenterCoarse[1] = self.oldcenterCoarse[1] + self.Coarseshift[1]
                if self.Coarseshift[1] < -1:
                    self.Coarseshift[1] = -1
                    self.newcenterCoarse[1] = self.oldcenterCoarse[1] + self.Coarseshift[1]
        else:
            print("DW not in!")

    def trackBubble(self):
        avgX = 0
        avgY = 0
        norm = 0
        mini = 0
        B = int(self.fineNodes[0] * 0.3)
        C = int(self.FinenodesPerCoarseCell[0])
        for i in range(self.fineNodes[0]):
            for j in range(self.fineNodes[1]):
                A = 0
                for k in range(self.fineNodes[2]):
                    #					norm += ((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1
                    #					avgX += i*(((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1)
                    #					avgY += j*(((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1)
                    #		self.newcenter = [int(avgX/norm), int(avgY/norm), self.oldcenter[2]]
                    A += self.fineSolver.state.M.get(i, j, k)[2] / (self.fineNodes[2] * self.mu)
                if A < mini:
                    mini = A
                    self.newcenter = [i, j, self.oldcenter[2]]
        shift = np.sqrt((self.newcenter[0] - self.oldcenter[0]) ** 2 + (self.newcenter[1] - self.oldcenter[1]) ** 2)
        if not (C < shift < B):
            self.newcenter = self.oldcenter
            self.Multiscaleloop()
        else:
            self.shift()



    def handle(self, state, tmpmesh1=None, tmpmesh2=None):
        if self.solver.state.t == 0:
            self.wholehandlestime = 0
            self.wholefineExch = 0
            self.wholecoarseExch = 0
        self.starthandle=time.time()
        #print("now in handle!")
        #print("simple external field value:"+str(self.fineSolver.state.fineExternalField.get(0,0,0)))
        #print("simple extrnal field uniform value:"+str(self.fineSolver.state.fineExternalField.uniform_value))
        # TODO: change that later, calculate it automatically!
        scale = self.mu / self.Ms
        while self.fineMagIsGiven == False:
            self.fillFineScale()
            #print("after filling finescale: "+str(self.fineSolver.state.M.get(2,2,0)))
            #print("after filling finescale coarse: "+str(self.solver.state.M.get(self.r_0[0])))
            #print("mu "+str(self.mu))
        #n= 1
        self.fineSolver.evolver._RungeKutta__controller.MAX_TIMESTEP = self.solver.state.h/4
        #self.solver.evolver._RungeKutta__controller.MAX_TIMESTEP = self.fineSolver.state.h*8
        growshrink = False

        if self.growshrink:
            if self.TrackingCond == 'SkyrmionEdgeT':
                if self.solver.state.shrinkGrowTr < 0:
                    self.shrinkright(self.solver.state.shrinkGrowTr )
                if self.solver.state.shrinkGrowTr > 0:
                    self.growright(self.solver.state.shrinkGrowTr )

                if self.solver.state.shrinkGrowTl > 0:
                    self.shrinkleft(self.solver.state.shrinkGrowTr )

                if self.solver.state.shrinkGrowTl < 0:
                    self.growleft(self.solver.state.shrinkGrowTr )

            if self.TrackingCond == 'SkyrmionEdgeB':
                if self.solver.state.shrinkGrowBr < 0:
                    self.shrinkright(self.solver.state.shrinkGrowBr )
                if self.solver.state.shrinkGrowBr > 0:
                    self.growright(self.solver.state.shrinkGrowBr )

                if self.solver.state.shrinkGrowBl > 0:
                    self.shrinkleft(self.solver.state.shrinkGrowBr )

                if self.solver.state.shrinkGrowBl < 0:
                    self.growleft(self.solver.state.shrinkGrowBr )
            if self.TrackingCond == 'SkyrmionEdgeR':
                if self.solver.state.shrinkGrowRb < 0:
                    self.growbottom(self.solver.state.shrinkGrowRb )
                if self.solver.state.shrinkGrowRb > 0:
                    self.shrinkbottom(self.solver.state.shrinkGrowRb )

                if self.solver.state.shrinkGrowRt > 0:
                    self.growtop(self.solver.state.shrinkGrowRt )

                if self.solver.state.shrinkGrowRt < 0:
                    self.shrinktop(self.solver.state.shrinkGrowRt )

            if self.TrackingCond == 'SkyrmionEdgeL':
                if self.solver.state.shrinkGrowLb < 0:
                    self.growbottom(self.solver.state.shrinkGrowLb )
                if self.solver.state.shrinkGrowLb > 0:
                    self.shrinkbottom(self.solver.state.shrinkGrowLb )
                if self.solver.state.shrinkGrowLt > 0:
                    self.growtop(self.solver.state.shrinkGrowLt )
                if self.solver.state.shrinkGrowLt < 0:
                    self.shrinktop(self.solver.state.shrinkGrowLt )

        if not self.t == 0:
            if self.Tracking == True:
                if self.TrackingCond == 'Vortex':
                    if self.solver.state.step % 10 == 0:
                        self.trackVortex()
                    else:
                    	self.Multiscaleloop()
            else:
                self.Multiscaleloop()
        else:
            self.Multiscaleloop()

        return 0 # self.fineSolver.state.M.to_numpy()

        #return 0

    def done(self):

        while self.fineMagIsGiven == False:
            self.fillFineScale()
        if not self.t == 0:
            #self.Tracking=True
            if self.Tracking == True:
                if self.TrackingCond == 'SkyrmionEdgeTR':
                    self.solver.state.TRinnerpoint = self.r_0
                    print("r_0:" +str(self.r_0))
                    print("TR innerpoint:" +str(self.solver.state.TRinnerpoint))

                if self.TrackingCond == 'SkyrmionEdgeBR':
                    self.solver.state.BRinnerpoint[1] = self.r_1[1]
                    print("r_1 : " + str(self.r_1))
                    print("BR innerpoint: " + str(self.r_1))

                if self.TrackingCond == 'SkyrmionEdgeTL':
                    self.solver.state.TLinnerpoint[0] = self.r_1[0]

                if self.TrackingCond == 'SkyrmionEdgeL':
                    self.trackSkyrmionEdgeL()
                    print("from L:" + str(self.Coarseshift))
                if self.TrackingCond == 'SkyrmionEdgeR':
                    self.trackSkyrmionEdgeR()
                    print("from R:" + str(self.Coarseshift))

                #    if self.solver.state.step % 10 == 0:
                if self.TrackingCond == 'SkyrmionEdgeT':
                    self.trackSkyrmionEdgeT()
                    print("from T:" + str(self.Coarseshift))
                if self.TrackingCond == 'SkyrmionEdgeB':
                    self.trackSkyrmionEdgeB()
                    print("from B:" + str(self.Coarseshift))
        #if not self.t == 0:

        self.stophandle=time.time()

        handletime = self.stophandle-self.starthandle
        self.wholehandlestime += handletime
        #print("handletime "+str(self.wholehandlestime))

        if not self.t == 0:
            self.handletime += (self.stophandle - self.starthandle)
            self.solvetime += (self.stopsolve - self.startsolve)
            #self.exchangetime += (self.stopexchange - self.startexchange)
            self.straytime += (self.stopfinestray - self.startfinestray)

            #print("whole: " + str(self.handletime))

            #print("solve: "+str(self.solvetime))
            #print("Exchange: "+str(self.exchangetime))
            #print("fineStray: "+str(self.straytime))
            #if (self.handletime>0):
            #    print("ratio fineStray whole: "+str(self.straytime/self.handletime))

        return 0
