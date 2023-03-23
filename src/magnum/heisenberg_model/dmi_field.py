# Copyright 2012, 2013 by the Micromagnum authors.
#
# This file is part of MicroMagnum.
#
# MicroMagnum is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MicroMagnum is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.

import magnum.module as module
import magnum.magneto as magneto
from magnum.mesh import VectorField, Field
from .constants import MU0

class FSDMIField(module.Module):
    def __init__(self):
        super(FSDMIField, self).__init__()

    def calculates(self):
        return ["H_dmi", "E_dmi"]

    def params(self):
        return ["Dx", "Dy", "Dz"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_fsdmi", 'EFFECTIVE_FIELD_ENERGY': "E_dmi"}

    def initialize(self, system):
        self.system = system
        Dx = VectorField(self.system.mesh); Dx.fill((0, 0, 0))
        Dy = VectorField(self.system.mesh); Dy.fill((0, 0, 0))
        Dz = VectorField(self.system.mesh); Dz.fill((0, 0, 0))
        setattr(self, 'Dx', Dx)
        setattr(self, 'Dy', Dy)
        setattr(self, 'Dz', Dz)

    def calculate(self, state, id):
        cache = state.cache

        if id == "H_dmi":
            if hasattr(cache, "H_fsdmi"): return cache.H_fsdmi
            H_fsdmi = cache.H_fsdmi = VectorField(self.system.mesh)
            H_fsdmi.fill((0,0,0))
            Dx = getattr(self, 'Dx')
            Dy = getattr(self, 'Dy')
            Dz = getattr(self, 'Dz')
            magneto.fs_dmi(self.system.mu, Dx, Dy, Dz, state.M, H_fsdmi)
            return H_fsdmi

        elif id == "E_dmi":
            return -MU0/2.0 * state.M.dotSum(state.H_fsdmi)

        else:
            raise KeyError("DMIField.calculate: Can't calculate %s", id)
