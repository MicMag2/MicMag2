#!/usr/bin/python
from magnum import *

world = World(
  RectangularMesh((100, 25, 2), (4e-9, 4e-9, 2e-9)),
  Body("bottom", Material.Py(), Cuboid((  0e-9, 0e-9, 0e-9), ( 400e-9, 100e-9, 2e-9))),
  Body("top", Material.Py(), Cuboid((0e-9, 0e-9, 2e-9), ( 400e-9, 100e-9, 4e-9))),
)

# Relax a domain wall
solver = create_solver(world, [StrayField, ExchangeField, InterlayerExchangeField], log=True)
solver.state.alpha = 0.5
solver.state["bottom"].M = (-8e5,   0, 0) # left side: left
solver.state["top"].M = (   0, 8e5, 0) # center: up
solver.relax()
mag = solver.state.M

# Apply current
solver = create_solver(world, [StrayField, ExchangeField,InterlayerExchangeField], log=True)
solver.state.M = mag
solver.j_offs = (4e11, 0, 0)
solver.addStepHandler(DataTableLog("domainwall.odt"), condition.EveryNthStep(20))
solver.addStepHandler(OOMMFStorage("domainwall", ["M", "H_stray"]), condition.EveryNthStep(20))

solver.solve(condition.Time(20e-9))
