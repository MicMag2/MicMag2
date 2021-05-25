from magnum.solver import StepHandler
import magnum.solver as solver
import magnum.magneto as magneto
import magnum.module as module
from constants import PI

class FSContinuityCheck(StepHandler):

	def __init__(self, treshold, M, Ms, nx, ny, nz):
		super(ContinuityCheck, self).__init__()
		self.treshold = treshold*PI
		self.system = module.system
		self.M = M
		self.Ms = Ms
		self.nx = nx
		self.ny = ny
		self.nz = nz
		print 'HM cc'

	def handle(self, state):
		self.r_0 = magneto.cosinecheck(self.nx, self.ny, self.nz, self.M, self.Ms, self.treshold).r_0
		self.r_1 = magneto.cosinecheck(self.nx, self.ny, self.nz, self.M, self.Ms, self.treshold).r_1
		print self.r_0
		print self.r_1
		return self.r_0, self.r_1


	def done(self):
		pass

