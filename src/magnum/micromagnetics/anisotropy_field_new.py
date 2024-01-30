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

class AnisotropyFieldNew(module.Module):
    def __init__(self):
        super(AnisotropyFieldNew, self).__init__()

    def calculates(self):
        return ["H_aniso_new", "E_aniso_new"]

    def params(self):
        return ["k_uniaxial1", "k_uniaxial2", "k_cubic1", "k_cubic2", "axisu", "axis1", "axis2"]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': "H_aniso_new", 'EFFECTIVE_FIELD_ENERGY': "E_aniso_new"}

    def initialize(self, system):
        self.system = system
        self.k_uniaxial1 = Field(self.system.mesh); self.k_uniaxial1.fill(0.0)
        self.k_uniaxial2 = Field(self.system.mesh); self.k_uniaxial2.fill(0.0)
        self.k_cubic1 = Field(self.system.mesh); self.k_cubic1.fill(0.0)
        self.k_cubic2 = Field(self.system.mesh); self.k_cubic2.fill(0.0)
        self.axisu = VectorField(self.system.mesh); self.axisu.fill((0.0, 0.0, 0.0))
        self.axis1 = VectorField(self.system.mesh); self.axis1.fill((0.0, 0.0, 0.0))
        self.axis2 = VectorField(self.system.mesh); self.axis2.fill((0.0, 0.0, 0.0))

    #def on_param_update(self, id):

        #if id in self.params() + ["Ms"]:
            #axis1, axis2 = self.axis1, self.axis2
            ##print axis1, axis2
            #k_uni, k_cub = self.k_uniaxial, self.k_cubic
            #Ms = self.system.Ms

            #def compute_none(state, H_aniso):
                #H_aniso.fill((0.0, 0.0, 0.0))
                #return 0.0

            #def compute_uniaxial(state, H_aniso):
                #return magneto.uniaxial_anisotropy(axis1, k_uni, Ms, state.M, H_aniso)

            #def compute_cubic(state, H_aniso):
                #return magneto.cubic_anisotropy(axis1, axis2, k_cub, Ms, state.M, H_aniso)

            #def compute_uniaxial_and_cubic(state, H_aniso):
                #tmp = VectorField(self.system.mesh)
                #E0 = magneto.uniaxial_anisotropy(axis1, k_uni, Ms, state.M, tmp)
                #E1 = magneto.cubic_anisotropy(axis1, axis2, k_cub, Ms, state.M, H_aniso)
                #state.cache.E_aniso_sum = E0 + E1
                #H_aniso.add(tmp)

            #fns = {(False, False): compute_none,
                   #( True, False): compute_uniaxial,
                   #(False,  True): compute_cubic,
                   #( True,  True): compute_uniaxial_and_cubic}

            #have_uni = not (k_uni.isUniform() and k_uni.uniform_value == 0.0)
            #have_cub = not (k_cub.isUniform() and k_cub.uniform_value == 0.0)
            #self.__compute_fn = fns[have_uni, have_cub]

    def calculate(self, state, id):
        if id == "H_aniso_new":
            if hasattr(state.cache, "H_aniso_new"): return state.cache.H_aniso_new
            H_aniso_new = state.cache.H_aniso_new = VectorField(self.system.mesh)
            state.cache.E_aniso_new_sum = magneto.anisotropy_new(self.axisu,
															self.k_uniaxial1,
															self.k_uniaxial2,
															self.axis1,
															self.axis2,
															self.k_cubic1,
															self.k_cubic2,
															self.system.Ms,
															state.M,
															H_aniso_new)
            return H_aniso_new

        elif id == "E_aniso_new":
            if not hasattr(state.cache, "E_aniso_new_sum"):
                self.calculate(state, "H_aniso_new")
            return state.cache.E_aniso_new_sum * self.system.mesh.cell_volume

        else:
            raise KeyError("AnisotropyFieldNew.calculate: Can't calculate %s", id)

