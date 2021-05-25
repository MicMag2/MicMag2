from __future__ import division
from magnum.solver import StepHandler
import math
import magnum.solver as solver
import magnum.magneto as magneto
import magnum.module as module
from magnum.mesh import RectangularMesh 
from magnum.micromagnetics.world import *
from constants import PI
from magnum.create_solver import create_solver
from magnum.heisenberg_model import *
from magnum.mesh import *
import magnum.solver.condition as condition
import numpy as np
from magnum.micromagnetics import StrayFieldCalculator

class ContinuityCheck(StepHandler):

	def __init__(self, coarseSolver, fineSolver, r_0 , r_1, FinenodesPerCoarseCell, folder, Tracking, TrackingCondition, fineMagIsGiven, DMI = False, mintrack=0.05, multiscaleLog = True):

		super(ContinuityCheck, self).__init__()
		self.DMI = DMI
		self.solver = coarseSolver
		self.world = coarseSolver.world
		self.mesh = coarseSolver.world.mesh
		self.NumNodes = self.mesh.getNumNodes()
		self.Delta = self.mesh.getDelta()
		self.mintrack = mintrack

		self.fineSolver = fineSolver
		self.fineWorld = fineSolver.world
		self.fineMesh = fineSolver.world.mesh
		self.fineNodes = self.fineMesh.getNumNodes()
		self.l = self.fineMesh.getDelta()

		self.r_0 = r_0
		self.r_1 = r_1

		self.body = self.world.bodies[0]
		self.Ms = self.body.material.Ms
		self.mu = self.body.material.mu
		if self.Ms == 0 or self.mu == 0:
			self.body = self.world.bodies[1]
			self.Ms = self.body.material.Ms
			self.mu = self.body.material.mu


		self.stray = StrayFieldCalculator(self.mesh)
		self.FinenodesPerCoarseCell = FinenodesPerCoarseCell

		self.zoomsize = [r_1[0] - r_0[0], r_1[1] - r_0[1], r_1[2] - r_0[2]] 
		self.oldcenter = [math.floor(self.fineNodes[0]/2) , math.floor(self.fineNodes[1]/2), math.floor(self.fineNodes[2]/2)]
		self.newcenter = [math.floor(self.fineNodes[0]/2) , math.floor(self.fineNodes[1]/2), math.floor(self.fineNodes[2]/2)]
		self.t = self.solver.state.t

		if self.NumNodes[2] == 1:
			if self.NumNodes[1] == 1:
				self.dim = 1
			else:
				self.dim = 2
		else:
			self.dim = 3
		self.Tracking = Tracking
		self.TrackingCond = TrackingCondition
		self.Fineshift = [0,0,0]
		self.fineMagIsGiven = fineMagIsGiven
		self.multiscaleLog = multiscaleLog
		self.folder = folder

	def oneDcenter(self, A, idx):
		return (A - self.r_0[idx])*self.FinenodesPerCoarseCell[idx] + int(self.FinenodesPerCoarseCell[idx]/2)

	def center(self, A):
		return [self.oneDcenter(A[0], 0), self.oneDcenter(A[1], 1), self.oneDcenter(A[2], 2)]

	def FinetoCoarse(self, q):
		Q = [0,0,0]
		for x in range(3):
			if q[x] >= 0:
				Q[x] = int(self.r_0[x] + math.floor(q[x]/self.FinenodesPerCoarseCell[x]))
			else:
				Q[x] = int(self.r_0[x] + math.ceil(q[x]/self.FinenodesPerCoarseCell[x]))
		return Q

	def chooseInterpolationCoordinates(self, i,j,k, I,J,K):
		xyz = [int(np.sign(i - self.oneDcenter(I,0))), int(np.sign(j - self.oneDcenter(J,1))), int(np.sign(k - self.oneDcenter(K,2)))]
		XYZ = [I, J, K]
		for x in range(3):
			if xyz[x] == 0:
				xyz[x] =  1
			if XYZ[x] == 0:
				xyz[x] =  1
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
		zero = [0,0,0]
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
		
	def linearInt(self, Q, fQ, x, axis = 0):
		A, B = self.center(Q[0]), self.center(Q[1])
		fA, fB = fQ[0], fQ[1]
		fx = [0,0,0]
		for coord in range(3):
			fx[coord] = fA[coord] + (fB[coord] - fA[coord])*(x[axis]-A[axis])/(B[axis]-A[axis])
		newQ = [0,0,0]
		for i in range(3):
			if i == axis:
				newQ[i] = self.FinetoCoarse(x)[i]
			else:
				newQ[i] = Q[0][i]
		return newQ, fx
		
	def bilinearInt(self, Q, fQ, x, axis0 = 0, axis1 = 1):
		newQ1, f1 = self.linearInt(Q[0], fQ[0], x, axis0)
		newQ2, f2 = self.linearInt(Q[1], fQ[1], x, axis0)
		newQ, fx = self.linearInt([newQ1, newQ2], [f1, f2], x, axis1)
		return newQ, fx
	
	def quadrilinearInt(self, Q, fQ, x, axis0 = 0, axis1 = 1, axis2 = 2):
		newQ1, f1 = self.bilinearInt(Q[0], fQ[0], x, axis0, axis1)
		newQ2, f2 = self.bilinearInt(Q[1], fQ[1], x, axis0, axis1)
		newQ, fx = self.linearInt([newQ1, newQ2], [f1, f2], x, axis2)
		return newQ, fx

	def Interpolate(self, i,j,k, I,J,K, Field, scale = 1):
		Q, fQ = self.chooseInterpolationCoordinates(i,j,k, I,J,K)
		if self.dim == 3:
			for l in range(2):
				for m in range(2):
					for n in range(2):
						Q_i = Q[l][m][n]
						fQ[l][m][n] = [(Field.get(*Q_i)[0])*scale, (Field.get(*Q_i)[1])*scale, (Field.get(*Q_i)[2])*scale] 
			newQ, fx = self.quadrilinearInt(Q, fQ, [i,j,k])
		elif self.dim == 2:
			for l in range(2):
				for m in range(2):
					Q_i = Q[l][m]
					fQ[l][m] = [(Field.get(*Q_i)[0])*scale, (Field.get(*Q_i)[1])*scale, (Field.get(*Q_i)[2])*scale]
			newQ, fx = self.bilinearInt(Q, fQ, [i,j,k])
		elif self.dim == 1:
			for l in range(2):
				Q_i = Q[l]
				fQ[l] = [(Field.get(*Q_i)[0])*scale, (Field.get(*Q_i)[1])*scale, (Field.get(*Q_i)[2])*scale] 
			newQ, fx = self.linearInt(Q, fQ, [i,j,k])
		return fx



	def mixedInterpolate(self, i,j,k, I,J,K, shift, fineshift, coarseField, fineField, scale = 1):
		Q, fQ = self.chooseInterpolationCoordinates(i,j,k, I,J,K)
		if self.dim == 3:
			for l in range(2):
				for m in range(2):
					for n in range(2):
						Q_i = Q[l][m][n]
						if (self.r_0[0] <= Q_i[0] + shift[0] <= self.r_1[0]) and (self.r_0[1] <=Q_i[1] + shift[1] <= self.r_1[1]):
							int_i, int_j, int_k = self.center(Q_i)
							fQ[l][m][n] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k )
						else:
							fQ[l][m][n] = [(coarseField.get(*Q_i)[0])*scale, (coarseField.get(*Q_i)[1])*scale, (coarseField.get(*Q_i)[2])*scale] 
			newQ, fx = self.quadrilinearInt(Q, fQ, [i,j,k])
		elif self.dim == 2:
			for l in range(2):
				for m in range(2):
					Q_i = Q[l][m]
					if (self.r_0[0] <=Q_i[0] + shift[0] <= self.r_1[0]) and (self.r_0[1] <=Q_i[1] + shift[1] <= self.r_1[1]):
						int_i, int_j, int_k = self.center(Q_i)
						fQ[l][m] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k )
					else:
						fQ[l][m] = [(coarseField.get(*Q_i)[0])*scale, (coarseField.get(*Q_i)[1])*scale, (coarseField.get(*Q_i)[2])*scale]
			newQ, fx = self.bilinearInt(Q, fQ, [i,j,k])
		elif self.dim == 1:
			for l in range(2):
				Q_i = Q[l]
				if (self.r_0[0] <= Q_i[0] + shift[0] <= self.r_1[0]) and (self.r_0[1] <= Q_i[1]+shift[1] <= self.r_1[1]):
					int_i, int_j, int_k = self.center(Q_i)
					fQ[l] = fineField.get(int_i + fineshift[0], int_j + fineshift[1], int_k )
				else:
					fQ[l] = [(coarseField.get(*Q_i)[0])*scale, (coarseField.get(*Q_i)[1])*scale, (coarseField.get(*Q_i)[2])*scale]
			newQ, fx = self.linearInt(Q, fQ, [i,j,k])
		return fx

	def fineStray(self):
		
		M_dummy = VectorField(self.mesh)
		H_Stray = VectorField(self.mesh)
		H_Stray.fill((0,0,0))
		for I in range (self.NumNodes[0]):
			for J in range (self.NumNodes[1]):
				for K in range (self.NumNodes[2]):
					if self.r_0[0] <= I <= self.r_1[0] and self.r_0[1] <= J <= self.r_1[1]and self.r_0[2] <= K <= self.r_1[2]:
						M_dummy.set(I,J,K,(0,0,0))
					else:
						M_dummy.set(I,J,K,(self.solver.state.M.get(I,J,K)))
		self.stray.calculate(M_dummy, H_Stray)
		for i in range(self.fineNodes[0]):
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.Interpolate(i,j,k, I,J,K, H_Stray)
					self.fineSolver.state.fineStrayField.set(i, j, k, (fx[0], fx[1], fx[2]))		
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

	def fineExch(self):
		left = 1
		right = 1
		top = 1
		bottom = 1
		below = 1
		above = 1
		if self.r_0[0] == 0:	left = 0 
		if self.r_1[0] == self.NumNodes[0] - 1:	right = 0 
		if self.r_0[1] == 0:	bottom = 0
		if self.r_1[1] == self.NumNodes[1] - 1:	top = 0 
		if self.r_0[2] == 0:	below = 0
		if self.r_1[2] == self.NumNodes[2] - 1:	above = 0
		self.dummymesh = RectangularMesh((self.fineNodes[0]+left+right, self.fineNodes[1]+top+bottom, self.fineNodes[2]+below+above), (self.l[0], self.l[1], self.l[2]))
		self.M_dummy = VectorField(self.dummymesh)
		self.M_dummy.fill((0,0,0))
		self.OtI_H_Exch = VectorField(self.dummymesh)
		self.OtI_H_Exch.fill((0,0,0))
		if self.DMI == True:
			self.OtI_H_DMI = VectorField(self.dummymesh)
			self.OtI_H_DMI.fill((0,0,0))
			dummyDx = VectorField(self.dummymesh)
			dummyDy = VectorField(self.dummymesh)
			dummyDz = VectorField(self.dummymesh)
			dummyDx.fill(self.fineSolver.state.Dx.get(0,0,0))
			dummyDy.fill(self.fineSolver.state.Dy.get(0,0,0))
			dummyDz.fill(self.fineSolver.state.Dz.get(0,0,0))
		self.mu_dummy = Field(self.dummymesh)
		self.mu_dummy.fill(self.body.material.mu)
		self.J_dummy = Field(self.dummymesh)
		self.J_dummy.fill(self.body.material.J)
		if left == 1:
			i = -1
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		if right == 1:
			i = self.fineNodes[0]
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		if bottom == 1:
			j = -1
			for i in range(self.fineNodes[0]):
				for k in range(self.fineNodes[2]):
					l,m,n =  i, j - bottom, k
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		if top == 1:
			j = self.fineNodes[1]
			for i in range(self.fineNodes[0]):
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		if below == 1:
			k = -1
			for i in range(self.fineNodes[0]):
				for j in range(self.fineNodes[1]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		if above == 1:
			k = self.fineNodes[2]
			for i in range(self.fineNodes[0]):
				for j in range(self.fineNodes[1]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.mixedInterpolate(i,j,k, I,J,K, [0,0,0], [0,0,0], self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
					self.M_dummy.set(i + left, j + bottom, k + below, (fx[0], fx[1], fx[2]))
		magneto.fs_exchange(self.mu_dummy, self.J_dummy, self.M_dummy, self.OtI_H_Exch)
		if self.DMI == True:
			magneto.fs_dmi(self.mu_dummy, dummyDx, dummyDy, dummyDz, self.M_dummy, self.OtI_H_DMI)
			

		if left == 1:
			i = 0
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))
						#print 'left: Exch=', self.fineSolver.state.fineExchField.get(i, j, k), 'DMI =', self.fineSolver.state.fineDMIField.get(i, j, k)
		if right == 1:
			i = self.fineNodes[0] - 1
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))
						#print 'right: Exch=', self.fineSolver.state.fineExchField.get(i, j, k), 'DMI =', self.fineSolver.state.fineDMIField.get(i, j, k)

		if bottom == 1:
			j = 0
			for i in range(self.fineNodes[0]):
				for k in range(self.fineNodes[2]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))
		if top == 1:
			j = self.fineNodes[1] - 1
			for i in range(self.fineNodes[0]):
				for k in range(self.fineNodes[2]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))
		if below == 1:
			k = 0
			for i in range(self.fineNodes[0]):
				for j in range(self.fineNodes[1]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))
		if above == 1:
			k = self.fineNodes[2] - 1
			for i in range(self.fineNodes[0]):
				for j in range(self.fineNodes[1]):
					self.fineSolver.state.fineExchField.set(i, j, k, self.OtI_H_Exch.get(i + left, j + bottom, k + below))
					if False:#self.DMI == True:
						self.fineSolver.state.fineDMIField.set(i, j, k, self.OtI_H_DMI.get(i + left, j + bottom, k + below))

	def cellAverage(self, field, scale = 1):
		cellVolume = self.FinenodesPerCoarseCell[0]*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[2]
		M = VectorField(self.mesh)
		M.fill((0,0,0))
		for I in range(self.r_0[0], self.r_1[0] + 1):
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					mag = [0,0,0]
					for l in range(self.FinenodesPerCoarseCell[0]):
						i = (I-self.r_0[0])*self.FinenodesPerCoarseCell[0] + l
						for m in range(self.FinenodesPerCoarseCell[1]):
							j = (J-self.r_0[1])*self.FinenodesPerCoarseCell[1] + m
							for n in range(self.FinenodesPerCoarseCell[2]):
								k = (K-self.r_0[2])*self.FinenodesPerCoarseCell[2] + n
								for x in range(3):
									mag[x] += (field.get(i,j,k)[x])*(scale/cellVolume)
					M.set(I,J,K, (mag[0], mag[1], mag[2]))
		return M
					
	def coarseStray(self):
		M = self.cellAverage(self.fineSolver.state.M, self.Ms/self.mu)
		self.stray.calculate(M, self.solver.state.coarseStrayField)
			
	def coarseExch(self):
		dx = 1	#(((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**-2) -1
		dy = 1	#(((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**-2) -1
		dz = 1	#(((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**-2) -1
		scale = self.body.material.mu/self.body.material.Ms
		H_tot = VectorField(self.mesh)
		H_tot.fill((0,0,0))
		H_totDMI = VectorField(self.mesh)
		H_totDMI.fill((0,0,0))

		if self.r_0[1] != 0:
			for I in range(self.r_0[0], self.r_1[0] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					wallmesh = RectangularMesh((self.zoomsize[0]+1, 2, self.zoomsize[2]+1), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					if self.DMI == True:
						H_DMI = VectorField(wallmesh)
						H_DMI.fill((0,0,0))
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for i in range(self.FinenodesPerCoarseCell[0]):
						for k in range(self.FinenodesPerCoarseCell[2]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0]*(I-self.r_0[0]) + i, 0, self.FinenodesPerCoarseCell[2]*(K-self.r_0[2]) + k)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[0]*self.FinenodesPerCoarseCell[2])#*((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**2)
					oldbottom = self.solver.state.M.get(I, self.r_0[1], K)
					M_dummy.set(I-self.r_0[0], 0, K-self.r_0[2], (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for I in range(self.r_0[0], self.r_1[0] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					H_tot.set(I, self.r_0[1]-1, K, H_Exch.get(I-self.r_0[0], 1, K-self.r_0[2]))
					if self.DMI == True:
						H_totDMI.set(I, self.r_0[1]-1, K, H_DMI.get(I-self.r_0[0], 1, K-self.r_0[2]))

		if self.r_1[1] != (self.NumNodes[1]-1):
			for I in range(self.r_0[0], self.r_1[0] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					wallmesh = RectangularMesh((self.zoomsize[0]+1, 2, self.zoomsize[2]+1), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					if self.DMI == True:
						H_DMI = VectorField(wallmesh)
						H_DMI.fill((0,0,0))
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for i in range(self.FinenodesPerCoarseCell[0]):
						for k in range(self.FinenodesPerCoarseCell[2]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0]*(I-self.r_0[0]) + i, self.fineNodes[1]-1, self.FinenodesPerCoarseCell[2]*(K-self.r_0[2]) + k)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[0]*self.FinenodesPerCoarseCell[2])#*((1/self.FinenodesPerCoarseCell[1]) + 1 /2)**2)
					oldbottom = self.solver.state.M.get(I, self.r_1[1], K)
					M_dummy.set(I-self.r_0[0], 0, K-self.r_0[2], (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for I in range(self.r_0[0], self.r_1[0] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					H_tot.set(I, self.r_1[1]+1, K, H_Exch.get(I-self.r_0[0], 1, K-self.r_0[2]))
					if self.DMI == True:
						H_totDMI.set(I, self.r_1[1]+1, K, H_DMI.get(I-self.r_0[0], 1, K-self.r_0[2]))

		if self.r_0[0] != 0:
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					wallmesh = RectangularMesh((2, self.zoomsize[1]+1, self.zoomsize[2]+1), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					if self.DMI == True:
						H_DMI = VectorField(wallmesh)
						H_DMI.fill((0,0,0))
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for j in range(self.FinenodesPerCoarseCell[1]):
						for k in range(self.FinenodesPerCoarseCell[2]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(0, self.FinenodesPerCoarseCell[1]*(J-self.r_0[1]) + j, self.FinenodesPerCoarseCell[2]*(K-self.r_0[2]) + k)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[2])#*((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**2)
					oldbottom = self.solver.state.M.get(self.r_0[0], J, K)
					M_dummy.set(0, J-self.r_0[1], K-self.r_0[2], (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, 1,1,1)
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					H_tot.set(self.r_0[0]-1, J, K, H_Exch.get(1, J-self.r_0[1], K-self.r_0[2]))
					if self.DMI == True:
						H_totDMI.set(self.r_0[0]-1, J, K, H_DMI.get(1, J-self.r_0[1], K-self.r_0[2]))

		if self.r_1[0] != (self.NumNodes[0]-1):
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					wallmesh = RectangularMesh((2, self.zoomsize[1]+1, self.zoomsize[2]+1), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					H_DMI = VectorField(wallmesh)
					H_DMI.fill((0,0,0))
					if self.DMI == True:
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for j in range(self.FinenodesPerCoarseCell[1]):
						for k in range(self.FinenodesPerCoarseCell[2]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(self.fineNodes[0]-1, self.FinenodesPerCoarseCell[1]*(J-self.r_0[1]) + j, self.FinenodesPerCoarseCell[2]*(K-self.r_0[2]) + k)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[2])#*((1/self.FinenodesPerCoarseCell[0]) + 1 /2)**2)  TEST
					oldbottom = self.solver.state.M.get(self.r_1[0], J, K)
					M_dummy.set(0, J-self.r_0[1], K-self.r_0[2], (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, 1,1,1) # TEST
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for K in range(self.r_0[2], self.r_1[2] + 1):
					H_tot.set(self.r_1[0]+1, J, K, H_Exch.get(1, J-self.r_0[1], K-self.r_0[2]))
					if self.DMI == True:
						H_totDMI.set(self.r_1[0]+1, J, K, H_DMI.get(1, J-self.r_0[1], K-self.r_0[2]))

		if self.r_0[2] != 0:
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for I in range(self.r_0[2], self.r_1[2] + 1):
					wallmesh = RectangularMesh((self.zoomsize[0]+1, self.zoomsize[1]+1, 2), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					if self.DMI == True:
						H_DMI = VectorField(wallmesh)
						H_DMI.fill((0,0,0))
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for i in range(self.FinenodesPerCoarseCell[0]):
						for j in range(self.FinenodesPerCoarseCell[1]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0]*(I-self.r_0[0]) + i, self.FinenodesPerCoarseCell[1]*(J-self.r_0[1]) + j, 0)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[0])#*((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**2)
					oldbottom = self.solver.state.M.get(I, J, self.r_0[2])
					M_dummy.set(I-self.r_0[0], J-self.r_0[1], 0, (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for I in range(self.r_0[0], self.r_1[0] + 1):
					H_tot.set(I, J, self.r_0[2]-1, H_Exch.get(I-self.r_0[0], J-self.r_0[1], 1))
					if self.DMI == True:
						H_totDMI.set(I, J, self.r_0[2]-1, H_DMI.get(I-self.r_0[0], J-self.r_0[1], 1))

		if self.r_1[2] != (self.NumNodes[2]-1):
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for I in range(self.r_0[0], self.r_1[0] + 1):
					wallmesh = RectangularMesh((self.zoomsize[0]+1, self.zoomsize[1]+1, 2), (self.Delta[0], self.Delta[1], self.Delta[2]))
					M_dummy = VectorField(wallmesh)
					M_dummy.fill((0,0,0))
					Ms_dummy = Field(wallmesh)
					Ms_dummy.fill(self.body.material.Ms)
					A_dummy = Field(wallmesh)
					A_dummy.fill(self.body.material.A)
					H_Exch = VectorField(wallmesh)
					H_Exch.fill((0,0,0))
					if self.DMI == True:
						H_DMI = VectorField(wallmesh)
						H_DMI.fill((0,0,0))
						wallDx = VectorField(wallmesh)
						wallDy = VectorField(wallmesh)
						wallDz = VectorField(wallmesh)
						wallDx.fill(self.solver.state.Dx.get(0,0,0))
						wallDy.fill(self.solver.state.Dy.get(0,0,0))
						wallDz.fill(self.solver.state.Dz.get(0,0,0))
					bottomside = [0,0,0]
					bottomrescaled = [0,0,0]
					for i in range(self.FinenodesPerCoarseCell[0]):
						for k in range(self.FinenodesPerCoarseCell[2]):
							for x in range(3):
								bottomside[x] += self.fineSolver.state.M.get(self.FinenodesPerCoarseCell[0]*(I-self.r_0[0]) + i, self.FinenodesPerCoarseCell[1]*(J-self.r_0[1]) + j, self.fineNodes[2]-1)[x]
					for x in range(3):		
						bottomrescaled[x] = bottomside[x]/(scale*self.FinenodesPerCoarseCell[1]*self.FinenodesPerCoarseCell[0])#*((1/self.FinenodesPerCoarseCell[2]) + 1 /2)**2)
					oldbottom = self.solver.state.M.get(I, J, self.r_1[2])
					M_dummy.set(I-self.r_0[0], J-self.r_0[1], 0, (bottomrescaled[0] - oldbottom[0], bottomrescaled[1] - oldbottom[1], bottomrescaled[2] - oldbottom[2]))
			magneto.exchange(Ms_dummy, A_dummy, M_dummy, H_Exch, dx, dy, dz)
			if self.DMI == True:
				magneto.dmi(Ms_dummy, wallDx, wallDy, wallDz, M_dummy, H_DMI)
			for J in range(self.r_0[1], self.r_1[1] + 1):
				for I in range(self.r_0[0], self.r_1[0] + 1):
					H_tot.set(I, J, self.r_1[2]+1, H_Exch.get(I-self.r_0[0], J-self.r_0[1], 1))
					if self.DMI == True:
						H_totDMI.set(I, J, self.r_1[2]+1, H_DMI.get(I-self.r_0[0], J-self.r_0[1], 1))
		self.solver.state.coarseExchField = H_tot
		if self.DMI == True:
			self.solver.state.coarseDMIField = H_totDMI

	def fillFineScale(self):
		for i in range(self.fineNodes[0]):
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					fx = self.Interpolate(i,j,k, I,J,K, self.solver.state.M, self.mu/self.Ms)
					self.fineSolver.state.M.set(i,j,k, (fx[0], fx[1], fx[2]))
		self.fineMagIsGiven = True
		
	def Multiscaleloop(self):
		self.h = self.solver.state.h
		self.t =self.solver.state.t
		self.fineSolver.state.fineStrayField.fill((0,0,0))
		self.fineSolver.state.fineExchField.fill((0,0,0))
		if self.DMI == True:
			self.fineSolver.state.fineDMIField.fill((0,0,0))		
		self.fineStray()
		self.fineExch()
		if not self.t == 0:
			if self.multiscaleLog:
				print '_____________________________'
				print 'fine'
			self.fineSolver.solve(condition.Time(self.t + (self.h)/2))
		else:
			self.f = open(self.folder + '/pos.txt', 'w')
			self.f.write(str([self.t, self.r_0, self.r_1]) + '\n')

		Mavg = self.cellAverage(self.fineSolver.state.M, self.Ms/self.mu)
		for I in range (self.NumNodes[0]):
			for J in range (self.NumNodes[1]):
				for K in range (self.NumNodes[2]):
					if not(self.r_0[0] <= I <= self.r_1[0] and self.r_0[1] <= J <= self.r_1[1]and self.r_0[2] <= K <= self.r_1[2]):
						Mavg.set(I,J,K,(self.solver.state.M.get(I,J,K)))
		self.solver.state.M = Mavg
		self.solver.state.coarseStrayField.fill((0,0,0))
		self.solver.state.coarseExchField.fill((0,0,0))
		if self.DMI == True:
			self.solver.state.coarseDMIField.fill((0,0,0))
		self.coarseExch()
		self.done()

	def shift(self):
		self.oldcenterCoarse = self.FinetoCoarse(self.oldcenter)
		self.newcenterCoarse = self.FinetoCoarse(self.newcenter)
		self.Coarseshift = [self.newcenterCoarse[0] - self.oldcenterCoarse[0], self.newcenterCoarse[1] - self.oldcenterCoarse[1], 0]
		self.newCoarseshift = [0,0,0]
		for x in range(2):
			if self.r_0[x] + self.Coarseshift[x] < 0:
				self.newCoarseshift[x] = - self.r_0[x]
			elif self.r_1[x] + self.Coarseshift[x] > self.NumNodes[x] - 1:
				self.newCoarseshift[x] = self.NumNodes[x] - 1 - self.r_1[x]
			else:
				self.newCoarseshift[x] = self.Coarseshift[x]
			self.r_0[x] += self.newCoarseshift[x]
			self.r_1[x] += self.newCoarseshift[x]
			self.Fineshift[x] = self.newCoarseshift[x]*self.FinenodesPerCoarseCell[x]
		new_M = VectorField(self.fineMesh)
		for i in range(self.fineNodes[0]):
			new_i = i + self.Fineshift[0]
			for j in range(self.fineNodes[1]):
				new_j = j + self.Fineshift[1]
				for k in range(self.fineNodes[2]):
					I, J, K = self.FinetoCoarse([i,j,k])
					if 0 <= new_i < self.fineNodes[0] and 0 <= new_j < self.fineNodes[1]:
						new_M.set(i,j,k, self.fineSolver.state.M.get(new_i, new_j, k))
					else:
						fx = self.mixedInterpolate(i,j,k, I,J,K, self.newCoarseshift, self.Fineshift, self.solver.state.M, self.fineSolver.state.M, scale = self.mu/self.Ms)
						new_M.set(i,j,k, (fx[0], fx[1], fx[2]))
		self.fineSolver.state.M = new_M
		#self.f.write(str([self.t, self.r_0, self.r_1]) + '\n')
		self.Multiscaleloop()

	def trackVortex(self):
		mini = 0
		maxi = 0
		coordmini = [0,0]
		coordmaxi = [0,0]
		B = int(self.fineNodes[0]*0.1)
		for i in range(int(self.oldcenter[0] - B), int(self.oldcenter[0] + B)):
			for j in range(int(self.oldcenter[1] - B), int(self.oldcenter[1] + B)):
				A = 0
				for k in range(self.fineNodes[2]):
					A += self.fineSolver.state.M.get(i,j,k)[2]/(self.fineNodes[2]*self.body.material.mu)
				if A < mini:
					mini = A
					coordmini = [i,j, self.oldcenter[2]]
				elif A > maxi:
					maxi = A
					coordmaxi = [i,j, self.oldcenter[2]]
		if (mini < -0.9) and (maxi > 0.9):
			x1 = min(coordmini[0], coordmaxi[0])
			x2 = max(coordmini[0], coordmaxi[0])
			y1 = min(coordmini[1], coordmaxi[1])
			y2 = max(coordmini[1], coordmaxi[1])
			self.newcenter = [x1 + int((x2 -x1)/2), y1 + int((y2 -y1)/2), self.newcenter[2]]			
		else:
			if maxi > abs(mini):
				self.newcenter = coordmaxi			
			elif abs(mini) > maxi:
				self.newcenter = coordmini			
		shift = np.sqrt((self.newcenter[0] - self.oldcenter[0])**2 +  (self.newcenter[1] - self.oldcenter[1])**2)
		C = int(self.FinenodesPerCoarseCell[0])
		if not(C < shift < B):
			self.newcenter =  self.oldcenter
			self.Multiscaleloop()
		else:
			self.shift()

	def trackDW(self):
		mini =  1
		B = int(self.fineNodes[0]*0.1)
		for i in range(int(self.oldcenter[0] - B), int(self.oldcenter[0] + B)):
			A = 0
			for j in range(self.fineNodes[1]):
				for k in range(self.fineNodes[2]):
					A += abs(self.fineSolver.state.M.get(i,j,k)[0])/(self.fineNodes[1]*self.fineNodes[2]*self.body.material.mu)
			if A < mini:
				mini = A
				self.newcenter[0] = i
		C = int(self.fineNodes[0]*self.mintrack)
		shift =  abs(self.newcenter[0] - self.oldcenter[0])
		if not(C < shift < B):
			self.newcenter[0] =  self.oldcenter[0]
			self.Multiscaleloop()
		else:
			if self.t == 0:
				self.Multiscaleloop()
			else:
				self.shift()



	def trackBubble(self):
		avgX = 0
		avgY = 0
		norm = 0
		mini = 0
		B = int(self.fineNodes[0]*0.3)
		C = int(self.FinenodesPerCoarseCell[0])
		for i in range(self.fineNodes[0]):
			for j in range(self.fineNodes[1]):
				A = 0
				for k in range(self.fineNodes[2]):
#					norm += ((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1
#					avgX += i*(((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1)
#					avgY += j*(((self.fineSolver.state.M.get(i,j,k)[2])/self.mu)-1)
#		self.newcenter = [int(avgX/norm), int(avgY/norm), self.oldcenter[2]]			
					A += self.fineSolver.state.M.get(i,j,k)[2]/(self.fineNodes[2]*self.mu)
				if A < mini:
					mini = A
					self.newcenter = [i,j, self.oldcenter[2]]
		shift = np.sqrt((self.newcenter[0] - self.oldcenter[0])**2 +  (self.newcenter[1] - self.oldcenter[1])**2)
		if not(C < shift < B):
			self.newcenter =  self.oldcenter
			self.Multiscaleloop()
		else:
			self.shift()


	def handle(self, state):
		while self.fineMagIsGiven == False:
			self.fillFineScale() 
		if not self.t == 0:
			if self.Tracking == True:
				if self.TrackingCond == 'Vortex':
					if self.solver.state.step % 10 == 0:
						self.trackVortex()
					else:
						self.Multiscaleloop()
				elif self.TrackingCond == 'DW':
					self.trackDW()
				elif self.TrackingCond == 'Bubble':
					if self.solver.state.step % 5 == 0:
						self.trackBubble()
					else:
						self.Multiscaleloop()
			else:
				self.Multiscaleloop()
		else:
			self.Multiscaleloop()

	
	def done(self):
		if self.multiscaleLog:
			print '_____________________________'
			print 'coarse'
		return 0
	
	
