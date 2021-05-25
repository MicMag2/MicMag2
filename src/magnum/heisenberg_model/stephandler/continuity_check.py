from magnum.solver import StepHandler
import magnum.solver as solver
import magnum.magneto as magneto
import magnum.module as module
from constants import PI

class FSContinuityCheck(StepHandler):

	def __init__(self, fineSolver, coarseSolver, r_0, r_1):
		super(FSContinuityCheck, self).__init__()
		self.fineSolver =  fineSolver
		self.coarseSolver = coarseSolver
		self.coarseWorld = coarseSolver.world
		self.coarseMesh = self.coarseWorld.mesh
		self.r_0 = r_0
		self.r_1 = r_1

	def handle(self, state):


	def done(self):
		pass

