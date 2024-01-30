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
        self.__func_DL = "MSTAngDep_a_DL_fn"
        self.__a_DL = "MSTAngDep_a_DL"
        setattr(self, self.__func_DL,  "1.")
        setattr(self, self.__a_DL,  None)
        self.__func_FL = "MSTAngDep_a_FL_fn"
        self.__a_FL = "MSTAngDep_a_FL"
        setattr(self, self.__func_FL,  "1.")
        setattr(self, self.__a_FL,  None)
        #raise NotImplementedError("The MacroSpinTorque module does not work yet.")

    def calculates(self):
        return ["dMdt_MSTAngleDep"]

    def params(self):
        return [self.__func_DL, self.__a_DL, self.__func_FL, self.__a_FL, "MSTAngDep_ref"]

    def properties(self):
        return {'LLGE_TERM': "dMdt_MSTAngleDep"}

    def initialize(self, system):
        self.system = system
        MSTAngDep_ref = (0,0,0)
        setattr(self, 'MSTAngDep_ref', MSTAngDep_ref)

    def calculate(self, state, id):
        cache = state.cache
        if id == "dMdt_MSTAngleDep":
            if hasattr(cache, "dMdt_MSTAngleDep"): return cache.dMdt_MSTAngleDep

            fn_DL          = getattr(self, self.__func_DL)
            a_DL           = getattr(self, self.__a_DL)
            fn_FL          = getattr(self, self.__func_FL)
            a_FL           = getattr(self, self.__a_FL)
            MSTAngDep_ref  = getattr(self, 'MSTAngDep_ref')           

            if not isinstance(MSTAngDep_ref, tuple):
                raise ValueError("MacroSpinTorqueAngleDep.calculates: MSTAngDep_ref is not a tuple.")

            if not isinstance(a_DL, (Field, float, int)): # with user function
                raise ValueError("MacroSpinTorqueAngleDep.calculates: No valid MSTAngDep_a_DL or MSTAngDep_a_DL_fn parameters set.")
            if not isinstance(a_FL, (Field, float, int)): # with user function
                raise ValueError("MacroSpinTorqueAngleDep.calculates: No valid MSTAngDep_a_FL or MSTAngDep_a_FL_fn parameters set.")

            dMdt_MSTAngleDep = cache.dMdt_MSTAngleDep = VectorField(self.system.mesh)
            dMdt_MSTAngleDep.fill((0.0,0.0,0.0))		#?????????????????????????????????

            # Calculate macro spin torque term due to Slonchewski
            nx, ny, nz = self.system.mesh.num_nodes
            dx, dy, dz = self.system.mesh.delta
            (rx,ry,rz) = MSTAngDep_ref

            if isinstance(a_DL, (float, int)):
                tmp = a_DL; a_DL = Field(self.system.mesh); a_DL.fill(tmp)
            if isinstance(a_FL, (float, int)):
                tmp = a_FL; a_FL = Field(self.system.mesh); a_FL.fill(tmp)

            magneto.fdm_slonchewski_angleDep(
              nx, ny, nz, dx, dy, dz, rx, ry, rz, state.t, #self.__do_precess,
              state.Ms, state.alpha, state.j, a_DL, a_FL, fn_DL, fn_FL,
              state.M, dMdt_MSTAngleDep
            )
            #print dMdt_MSTAngleDep.get(100,100,0)
            return dMdt_MSTAngleDep

        else:
            raise KeyError("MacroSpinTorqueAngleDep.calculate: Can't calculate %s", id)
