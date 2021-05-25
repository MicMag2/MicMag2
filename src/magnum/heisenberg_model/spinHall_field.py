import magnum.module as module
import magnum.magneto as magneto
from magnum.mesh import VectorField, Field
from .constants import *

class FSSpinHallField(module.Module):
    def __init__(self):
        super(FSSpinHallField, self).__init__()

    def calculates(self):
        return ["H_SHE", "E_SHE"]

    def params(self):
        return ["Theta", "thickness", "Ms"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_exch", 'EFFECTIVE_FIELD_ENERGY': "E_exch"}

    def initialize(self, system):
        self.system = system
        self.thickness = self.system.thickness
        self.Theta = self.system.Theta

    def calculate(self, state, id):
        cache = state.cache

        if id == "H_SHE":
            if hasattr(cache, "H_SHE"): return cache.H_SHE
            H_SHE = cache.H_SHE = VectorField(self.system.mesh)
            self.mu = Field(self.system.mesh); self.mu.fill(self.system.mu)
            self.J = VectorField(self.system.mesh); self.J.fill(state.pj[0])
            magneto.fs_spinhall(self.mu, self.J, self.thickness, self.Ms, state.M, self.Theta, H_SHE)
            return H_SHE

        elif id == "E_SHE":
            return -MU0/2.0 * self.system.mesh.cell_volume * state.M.dotSum(state.H_SHE)

        else:
            raise KeyError("FSSpinHallField.calculate: Can't calculate %s", id)
