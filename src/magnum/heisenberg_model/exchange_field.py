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

class FSExchangeField(module.Module):
    def __init__(self):
        super(FSExchangeField, self).__init__()

    def calculates(self):
        return ["H_exch", "E_exch"]

    def params(self):
        return ["J", "mu"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_exch", 'EFFECTIVE_FIELD_ENERGY': "E_exch"}

    def initialize(self, system):
        self.system = system
        self.J = Field(self.system.mesh); self.J.fill(0.0)
        self.mu = Field(self.system.mesh); self.mu.fill(0.0)
        self.__peri_x = system.mesh.periodic_bc[0].find("x") != -1
        self.__peri_y = system.mesh.periodic_bc[0].find("y") != -1
        self.__peri_z = system.mesh.periodic_bc[0].find("z") != -1

    def calculate(self, state, id):
        cache = state.cache

        if id == "H_exch":
            if hasattr(cache, "H_exch"): return cache.H_exch
            H_exch = cache.H_exch = VectorField(self.system.mesh)
            H_exch.fill((0,0,0))
            #H_exch.fill((0,0,0))
            #nx, ny, nz = self.system.mesh.num_nodes
            #dx, dy, dz = self.system.mesh.delta
            #bcx, bcy, bcz = self.__peri_x, self.__peri_y, self.__peri_z
#            self.mu = Field(self.system.mesh); self.mu.fill(self.system.mu)     #TODO Fix this for multilayer samples
#            self.J = Field(self.system.mesh); self.J.fill(self.system.J)
            magneto.fs_exchange(self.system.mu, self.system.J, state.M, H_exch)
#            print('_________________________')
#            print(H_exch.get(0, 0, 0))
            #print "H_exch ", H_exch.average()
            #print "J ", self.system.J.average()
            #print "self.system.mu " , self.system.mu.average()
            #print "state.M.avg ", state.M.average()
            return H_exch

        elif id == "E_exch":
            return -MU0/2.0 * state.M.dotSum(state.H_exch)

        else:
            raise KeyError("FSExchangeField.calculate: Can't calculate %s", id)
