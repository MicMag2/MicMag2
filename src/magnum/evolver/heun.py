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

from magnum.mesh import VectorField

from .evolver import Evolver

class Heun(Evolver):
    def __init__(self, mesh, step_size):
        super(Heun, self).__init__(mesh)
        self.step_size = float(step_size)

    def evolve(self, state, t_max):
        #h = self.step_size
        #s0 = state
        dydt1 = state.differentiate()#; s0.flush_cache()
        #stmp = s0.clone(); stmp.substep=1; stmp.t = state.t

        state.y.add(dydt1, self.step_size)#; s0.t +=h
        #stmp = state.clone() ;
        #stmp.y.add(dydt1, self.step_size)

        state.finish_step()
        dydt2 = state.differentiate2()#; s0.flush_cache()
        #dydt2 = state.differentiate2()#; s0.flush_cache()
        #print("now diff2")
        state.y.add(dydt1, -self.step_size/2)
        state.y.add(dydt2, self.step_size/2)

        state.t += self.step_size
        state.h = self.step_size
        state.step += 1
        state.substep = 0
        state.flush_cache()
        state.finish_step()
        return state
