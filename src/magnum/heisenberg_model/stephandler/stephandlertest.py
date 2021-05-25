from magnum.solver import StepHandler



class test(StepHandler):

	def __init__(self):
		super(test, self).__init__()

	def handle(self, state):
		print 'funziono'

	def done(self):
		pass

