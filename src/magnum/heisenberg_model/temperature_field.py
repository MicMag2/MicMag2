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
#import magnum.evolver as evolver
from magnum.logger import logger
from .constants import MU0, BOLTZMANN_CONSTANT, GYROMAGNETIC_RATIO
import numpy as np
import time

from .io import writeOMF

class FSTemperatureField(module.Module):
    def __init__(self):
        super(FSTemperatureField, self).__init__()

    def calculates(self):
        return ["H_th", "E_th"]

    def params(self):
        return ["kelv_seed"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_th", 'EFFECTIVE_FIELD_ENERGY': "E_th"}

    def initialize(self, system):
        self.system = system
        self.kelv_seed = None
        #if not isinstance(self.evolver, evolver.Euler): raise NotImplementedError("TemperatureField: Only implemented for Euler solver!")
        logger.warn("The Temperature Field only works properly for fixed step sizes. This is currently only fulfilled by the Euler or Heun evolver! Furthermore, CPU and GPU currently use different random generators!")

    def calculate(self, state, id):
        cache = state.cache

        # if no seed given, generate seed from system time
        if state.kelv_seed == None:
            state.kelv_seed = time.time()

        if id == "H_th":
            if hasattr(cache, "H_th"): return cache.H_th
            H_th = cache.H_th = VectorField(self.system.mesh)
            H_th.fill((0,0,0))

            nx, ny, nz = self.system.mesh.num_nodes
            dx, dy, dz = self.system.mesh.delta
            #bcx, bcy, bcz = self.__peri_x, self.__peri_y, self.__peri_z

            magneto.fdm_temperature(nx,ny,nz,dx,dy,dz,self.system.Ms,self.system.alpha,state.kelv,state.h,state.step,state.kelv_seed,H_th)
	    #$print(H_th.get(40, 40, 0))
            return H_th

        elif id == "E_th":
            return -MU0/2.0 * self.system.mesh.cell_volume * state.M.dotSum(state.H_th)

        else:
            raise KeyError("TemperatureField.calculate: Can't calculate %s", id)


