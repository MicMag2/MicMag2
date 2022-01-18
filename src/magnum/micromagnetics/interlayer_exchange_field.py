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
from magnum.mesh import VectorField, Field, RectangularMesh
from .constants import MU0

class InterlayerExchangeField(module.Module):
    def __init__(self):
        super(InterlayerExchangeField, self).__init__()

    def calculates(self):
        return ["H_intexch", "E_intexch"]

    def params(self):
        return ["intExchMatrix"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_intexch", 'EFFECTIVE_FIELD_ENERGY': "E_intexch"}

    def initialize(self, system):
        self.system = system
        self.__peri_x = system.mesh.periodic_bc[0].find("x") != -1
        self.__peri_y = system.mesh.periodic_bc[0].find("y") != -1
        self.__peri_z = system.mesh.periodic_bc[0].find("z") != -1
        intExchMatrix = VectorField(self.system.mesh); intExchMatrix.fill((0, 0, 0))
        setattr(self, 'intExchMatrix', intExchMatrix)

    def calculate(self, state, id):
        cache = state.cache

        if id == "H_intexch":
            if hasattr(cache, "H_intexch"): return cache.H_intexch
            H_intexch = cache.H_intexch = VectorField(self.system.mesh)
            H_intexch.fill((0,0,0))

            #nx, ny, nz = self.system.mesh.num_nodes
            #dx, dy, dz = self.system.mesh.delta
            #bcx, bcy, bcz = self.__peri_x, self.__peri_y, self.__peri_z

            #setting up parameters
            intExchMatrix = getattr(self, 'intExchMatrix')

            magneto.interlayerExchange(self.system.Ms, self.intExchMatrix, state.M, H_intexch)
            #print H_intexch.get(100,25,0)
            return H_intexch

        elif id == "E_intexch":
            return -MU0/2.0 * self.system.mesh.cell_volume * state.M.dotSum(state.H_intexch)

        else:
            raise KeyError("InterlayerExchangeField.calculate: Can't calculate %s", id)
