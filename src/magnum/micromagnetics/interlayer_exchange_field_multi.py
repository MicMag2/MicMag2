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

class MultiInterlayerExchangeField(module.Module):
    def __init__(self):
        super(MultiInterlayerExchangeField, self).__init__()

    def calculates(self):
        return ["H_intexch_multi", "E_intexch_multi"]

    def params(self):
        return ["intExchPat"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_intexch_multi", 'EFFECTIVE_FIELD_ENERGY': "E_intexch_multi"}

    def initialize(self, system):
        self.system = system
        self.intExchPat =  []
        self.__peri_x = system.mesh.periodic_bc[0].find("x") != -1
        self.__peri_y = system.mesh.periodic_bc[0].find("y") != -1
        self.__peri_z = system.mesh.periodic_bc[0].find("z") != -1
	self.patternMatrix = VectorField(RectangularMesh((len(self.system.intExchPat), 1, 1), (1e-9, 1e-9, 1e-9)))

    def calculate(self, state, id):
        cache = state.cache

        if id == "H_intexch_multi":
            if hasattr(cache, "H_intexch_multi"): return cache.H_intexch_multi
            H_intexch_multi = cache.H_intexch_multi = VectorField(self.system.mesh)
            H_intexch_multi.fill((0,0,0))

            #nx, ny, nz = self.system.mesh.num_nodes
            #dx, dy, dz = self.system.mesh.delta
            #bcx, bcy, bcz = self.__peri_x, self.__peri_y, self.__peri_z

            #setting up parameters
            (nx, ny, nz) = self.patternMatrix.mesh.getNumNodes()

            if not nx == len(self.system.intExchPat):

                self.patternMatrix = VectorField(RectangularMesh((len(self.system.intExchPat), 1, 1), (1e-9, 1e-9, 1e-9)))

                for entryNum1 in range(len(self.system.intExchPat)):
                    for entryNum2 in range(len(self.system.intExchPat)):

                        if entryNum1 == entryNum2:
                            if self.system.intExchPat[entryNum1][0] == self.system.intExchPat[entryNum1][1]:
                                raise ValueError("InterlayerExchangeField_multi.calculate: pattern contains double entries.")
                            continue

                        if ((self.system.intExchPat[entryNum1][0] == self.system.intExchPat[entryNum2][0]) and (self.system.intExchPat[entryNum1][1] == self.system.intExchPat[entryNum2][1])) or ((self.system.intExchPat[entryNum1][1] == self.system.intExchPat[entryNum2][0]) and (self.system.intExchPat[entryNum1][0] == self.system.intExchPat[entryNum2][1])):
                            raise ValueError("InterlayerExchangeField_multi.calculate: pattern contains double entries.")
     
                for entryNum in range(len(self.system.intExchPat)):
                    self.patternMatrix.set(entryNum, 0, 0, (self.system.intExchPat[entryNum][0], self.system.intExchPat[entryNum][1], self.system.intExchPat[entryNum][2]))                

            magneto.interlayerExchange_multi(self.system.Ms, self.patternMatrix, state.M, H_intexch_multi)
            #print H_intexch_multi.get(100,25,0)
            return H_intexch_multi

        elif id == "E_intexch_multi":
            return -MU0/2.0 * self.system.mesh.cell_volume * state.M.dotSum(state.H_intexch_multi)

        else:
            raise KeyError("InterlayerExchangeField_multi.calculate: Can't calculate %s", id)

