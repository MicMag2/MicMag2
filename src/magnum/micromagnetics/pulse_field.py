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
from magnum.mesh import VectorField, Field
from magnum.logger import logger
from math import sin
import numpy as np

# PulseField
class PulseField(module.Module):
    def __init__(self, var_id):
        super(PulseField, self).__init__()

        # Parameters
        self.__offs = var_id + "_offs"
        self.__amp  = var_id + "_amp"
        self.__sigma = var_id + "_sigma"
        self.__mpv = var_id + "_mpv"
        self.__func = var_id + "_fn"

        setattr(self, self.__offs,  (0.0, 0.0, 0.0))
        setattr(self, self.__amp,   (0.0, 0.0, 0.0))
        setattr(self, self.__sigma,  (0.0, 0.0, 0.0))
        setattr(self, self.__mpv, (0.0, 0.0, 0.0))
        setattr(self, self.__func,  None)

        # Generated model variables
        self.__var = var_id

    def calculates(self):
        return [self.__var]

    def params(self):
        return [self.__offs, self.__amp, self.__sigma, self.__mpv, self.__func]

    def properties(self):
        return {'EFFECTIVE_FIELD_TERM': self.__var}
        
    def initialize(self, system):
        self.system = system
        logger.info("%s: Providing model variable %s, parameters are %s" % (self.name(), self.__var, ", ".join(self.params())))

    def calculate(self, state, id):
        if id == self.__var:
            # Get parameters...
            offs  = getattr(self, self.__offs)
            amp   = getattr(self, self.__amp)
            sigma  = getattr(self, self.__sigma)
            mpv = getattr(self, self.__mpv)
            fn    = getattr(self, self.__func)

            # Calculate field 'A'.
            t = state.t
            if fn: # with user function
                if any(x != (0.0, 0.0, 0.0) for x in (amp, sigma, mpv, offs)):
                    raise ValueError("AlternatingField.calculates: If %s is defined, the parameters %s, %s, %s and %s must be zero vectors, i.e. (0.0, 0.0, 0.0)" % (self.__func, self.__offs, self.__amp, self.__sigma, self.__mpv))
                # call user function
                A = fn(t)
            else:
                # with 'offs', 'amp', 'freq', 'phase' parameters
                A = (offs[0] + amp[0] * np.exp(-((t - mpv[0])**2)/(2*sigma[0]*sigma[0])),
                     offs[1] + amp[1] * np.exp(-((t - mpv[1])**2)/(2*sigma[1]*sigma[1])),
                     offs[2] + amp[2] * np.exp(-((t - mpv[2])**2)/(2*sigma[2]*sigma[2])),)

            # Convert 3-vector to VectorField if necessary.
            if isinstance(A, tuple):
                tmp = A; A = VectorField(self.system.mesh); A.fill(tmp)

            # Return field 'A'
            return A

        else:
            raise KeyError("AlternatingField.calculates: Can't calculate %s", id)
