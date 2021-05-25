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
from .constants import MU0, H_BAR, ELECTRON_CHARGE

# void fdm_slonchewski(
#       int dim_x, int dim_y, int dim_z,
#       double delta_x, double delta_y, double delta_z,
#       double a_j,
#       const VectorMatrix &p, // spin polarization
#       const Matrix &Ms,
#       const Matrix &alpha,
#       const VectorMatrix &M,
#       VectorMatrix &dM
# );

class FSMacroSpinTorque(module.Module):
    def __init__(self, do_precess = True):
        super(FSMacroSpinTorque, self).__init__()
        self.__do_precess = do_precess
        #raise NotImplementedError("The MacroSpinTorque module does not work yet.")

    def calculates(self):
        return ["dMdt_MST"]

    def params(self):
        return ["MST_xi", "MST_alpha_hall"]

    def properties(self):
        return {'LLGE_TERM': "dMdt_MST"}

    def initialize(self, system):
        self.system = system
        MST_xi = Field(self.system.mesh); MST_xi.fill(0.0)
        MST_alpha_hall = Field(self.system.mesh); MST_alpha_hall.fill(0.0)
        setattr(self, 'MST_xi', MST_xi)
        setattr(self, 'MST_alpha_hall', MST_alpha_hall)

    def calculate(self, state, id):
        cache = state.cache
        if id == "dMdt_MST":
            if hasattr(cache, "dMdt_MST"): return cache.dMdt_MST

            MST_xi             = getattr(self, 'MST_xi')
            MST_alpha_hall     = getattr(self, 'MST_alpha_hall')
            t = state.t

            dMdt_MST = cache.dMdt_MST = VectorField(self.system.mesh)
            dMdt_MST.fill((0.0,0.0,0.0))		#?????????????????????????????????
            #print "in", dMdt_MST.get(100,100,0)

            # Calculate macro spin torque term due to Slonchewski
            nx, ny, nz = self.system.mesh.num_nodes
            dx, dy, dz = self.system.mesh.delta

            if isinstance(MST_xi, (float, int)):
                tmp = MST_xi; MST_xi = Field(self.system.mesh); MST_xi.fill(tmp)

            if isinstance(MST_alpha_hall, (float, int)):
                tmp = MST_alpha_hall; MST_alpha_hall = Field(self.system.mesh); MST_alpha_hall.fill(tmp)

            self.Ms = Field(state.M.mesh)
            self.Ms.fill((state.mu.get(0,0,0)/((state.M.mesh.getDelta()[0])**3)))            
 
            magneto.FSfdm_slonchewski(
              nx, ny, nz, dx, dy, dz, #self.__do_precess,
              self.Ms, state.j, state.alpha, MST_xi, MST_alpha_hall,
              state.M, dMdt_MST, state.mu
            )
            #print dMdt_MST.get(1,1,0)
            return dMdt_MST

        else:
            raise KeyError("MacroSpinTorque.calculate: Can't calculate %s", id)
