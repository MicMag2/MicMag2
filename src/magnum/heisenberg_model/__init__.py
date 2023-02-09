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

## module system and solver for micromagnetics
from .heisenberg_model import HeisenbergModel
from .heisenberg_model_solver import HeisenbergModelSolver
from magnum.create_solver import create_solver

## constants and modules
from .constants import MU0, H_BAR, ELECTRON_CHARGE, MU_BOHR, GYROMAGNETIC_RATIO, BOLTZMANN_CONSTANT
from .landau_lifshitz_gilbert import LandauLifshitzGilbert
from .exchange_field import FSExchangeField
from .stray_field import FSStrayField, FSStrayFieldCalculator, FSStrayFieldTensor
from .external_field import ExternalField
from .anisotropy_field import FSAnisotropyField
from .homogeneous_field import HomogeneousField, HomogeneousCurrent
from .spin_torque import SpinTorque
from .alternating_field import AlternatingField
from .alternating_current import AlternatingCurrent
from .simple_field import SimpleExternalField, SimpleVectorField
from .pulse_current import PulseCurrent
from .spinHall_field import FSSpinHallField
from .pulse_field import PulseField
from .macro_spin_torque import FSMacroSpinTorque
from .dmi_field import FSDMIField
from .temperature_field import TemperatureField

__all__ = [
    "HeisenbergModel", "HeisenbergModelSolver", "create_solver",
    "MU0", "H_BAR", "ELECTRON_CHARGE", "MU_BOHR", "GYROMAGNETIC_RATIO", "BOLTZMANN_CONSTANT",
    "LandauLifshitzGilbert", "FSExchangeField", "FSStrayField", "FSStrayFieldCalculator","FSStrayFieldTensor",
    "ExternalField", "FSAnisotropyField", "HomogeneousField", "HomogeneousCurrent",
    "SpinTorque", "AlternatingField", "AlternatingCurrent", "SimpleExternalField", "SimpleVectorField", "FSSpinHallField", "PulseCurrent", "PulseField", "FSMacroSpinTorque", "FSDMIField", "TemperatureField"
]

## submodules
from . import io
from . import world
from . import stephandler
from . import toolbox

__all__.extend(world.__all__ + stephandler.__all__ + toolbox.__all__ + io.__all__)
from .world import *
from .stephandler import *
from .toolbox import *
from .io import *
