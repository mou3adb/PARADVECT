"""
To have the full formulation of e(t) = ebar + ehat*sin(2pi t/T + phi), we need
to compute e(t) at three values, 0, T/2 and T/4. For this, we write this
small script below.
"""
import time
import datetime
import numpy as np

from flow import Flow

from particle import Particle

from automatic_search_scripts import automatic_search

from results_fsi.results_fsi import *
#==============================================================================
printIt = False

parent_folder = r'../'

case = 'fixed_cylinder'

flow_path = parent_folder + 'flows/%s_atRe%d' % (case,Re)

elements_path   = parent_folder + 'elements/%s'   % case
geometries_path = parent_folder + 'geometries/%s' % case

flow = Flow()
flow.load_flow(flow_path, elements_path, geometries_path)

diameter = 0.1
density  = 2.

birth = 0
lifetime = 200

pos0 = np.array([-1.99, 0.01])
u0 = None

factor = 1

# Ready to start

increment = 0.01
precision = diameter*0.001

p = Particle(diameter, density, birth, lifetime, pos0, u0)
p.compute_trajectory(flow, factor, printIt)


e, limits = automatic_search(p, increment, precision,
                             flow, factor, printIt)

