from .alternating_field import FSAlternatingField

class PulseCurrent(FSAlternatingField):
	def __init__(self, var_id = "pj"):
		super(PulseCurrent, self).__init__(var_id)

		self.__amp = var_id + "_amp"
		self.__time = var_id + "_time"

		setattr(self, self.__time, 0)
		setattr(self, self.__amp,   (0.0, 0.0, 0.0))

	def calculates(self):
		return [self.__var]

	def params(self):
		return [self.__amp, self.__time]

	def properties(self):
		return {}

	def initialize(self, system):
		self.system = system
		logger.info("%s: Providing model variable %s, parameters are %s" % (self.name(), self.__var, ", ".join(self.params())))

	def calculate(self, state, id):
		if id == self.__var:
			# Get parameters...
			amp   = getattr(self, self.__amp)
			time    = getattr(self, self.__time)

			# Calculate field 'A'.
			t = state.t
			if t < self.time:
				A = (amp[0], amp[1], amp[2])
			else:
				A = (0, 0, 0)

			# Convert 3-vector to VectorField if necessary.
			if isinstance(A, tuple):
				tmp = A; A = VectorField(self.system.mesh); A.fill(tmp)

				# Return field 'A'
				return A

			else:
				raise KeyError("PulseCurrent.calculates: Can't calculate %s", id)

