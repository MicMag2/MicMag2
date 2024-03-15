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

class MacroSpinTorqueAngleDep(module.Module):
    def __init__(self, do_precess = True):
        super(MacroSpinTorqueAngleDep, self).__init__()
        self.__do_precess = do_precess
        #raise NotImplementedError("The MacroSpinTorque module does not work yet.")

    def calculates(self):
        return ["dMdt_MSTAngleDep"]

    def params(self):
        return ["MSTAngleDep_xi", "MSTAngleDep_alpha_hall", "MSTAngleDep_a_DL_fn", "MSTAngleDep_a_FL_fn"]

    def properties(self):
        return {'LLGE_TERM': "dMdt_MSTAngleDep"}

    def initialize(self, system):
        self.system = system
        MSTAngleDep_xi = Field(self.system.mesh); MSTAngleDep_xi.fill(0.0)
        MSTAngleDep_alpha_hall = Field(self.system.mesh); MSTAngleDep_alpha_hall.fill(0.0)
        
        MSTAngleDep_a_DL_fn = "1."
        MSTAngleDep_a_FL_fn = "1."
        
        setattr(self, 'MSTAngleDep_a_DL_fn', MSTAngleDep_a_DL_fn)
        setattr(self, 'MSTAngleDep_a_FL_fn', MSTAngleDep_a_FL_fn)        
        setattr(self, 'MSTAngleDep_xi', MSTAngleDep_xi)
        setattr(self, 'MSTAngleDep_alpha_hall', MSTAngleDep_alpha_hall) 
        

    def calculate(self, state, id):
        cache = state.cache
        if id == "dMdt_MSTAngleDep":
            if hasattr(cache, "dMdt_MSTAngleDep"): return cache.dMdt_MSTAngleDep

            MSTAngleDep_a_DL_fn        = getattr(self, "MSTAngleDep_a_DL_fn")
            MSTAngleDep_a_FL_fn        = getattr(self, "MSTAngleDep_a_FL_fn")       
            
            MSTAngleDep_xi             = getattr(self, 'MSTAngleDep_xi')
            MSTAngleDep_alpha_hall     = getattr(self, 'MSTAngleDep_alpha_hall')
            t = state.t
            
            
            #print("xi", MSTAngleDep_xi.get(0,0,0))
            #print("al", MSTAngleDep_alpha_hall.get(0,0,0))
            #print("DL", MSTAngleDep_a_DL_fn)
            #print("FL", MSTAngleDep_a_FL_fn)
            

            dMdt_MSTAngleDep = cache.dMdt_MSTAngleDep = VectorField(self.system.mesh)
            dMdt_MSTAngleDep.fill((0.0,0.0,0.0))		#?????????????????????????????????

            # Calculate macro spin torque term due to Slonchewski
            nx, ny, nz = self.system.mesh.num_nodes
            dx, dy, dz = self.system.mesh.delta
            
            if isinstance(MSTAngleDep_xi, (float, int)):
                tmp = MSTAngleDep_xi; MSTAngleDep_xi = Field(self.system.mesh); MSTAngleDep_xi.fill(tmp)

            if isinstance(MSTAngleDep_alpha_hall, (float, int)):
                tmp = MSTAngleDep_alpha_hall; MSTAngleDep_alpha_hall = Field(self.system.mesh); MSTAngleDep_alpha_hall.fill(tmp)            

            magneto.fdm_slonchewski_angleDep(
              nx, ny, nz, dx, dy, dz, state.t, #self.__do_precess,
              state.Ms, state.alpha, state.j, MSTAngleDep_xi, MSTAngleDep_alpha_hall, MSTAngleDep_a_DL_fn, MSTAngleDep_a_FL_fn, state.M, dMdt_MSTAngleDep
            )
            #print dMdt_MSTAngleDep.get(100,100,0)
            return dMdt_MSTAngleDep

        else:
            raise KeyError("MacroSpinTorqueAngleDep.calculate: Can't calculate %s", id)
