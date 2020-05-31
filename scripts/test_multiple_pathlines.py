"""
This is a script to test a particle simulation.

Please change any ..._path to your corresponding file path.

"""
import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as pp

from flow import Flow

from test_multiple_pathlines_scripts import prepare_initial_positions, spread_pathlines

from animationpathlines import AnimationPathlines
#==============================================================================
# Prepare flow
parent_folder = r'../'

case = 'potential'

flow_path       = parent_folder + 'flows/%s' % case
elements_path   = parent_folder + 'elements/%s' % case
geometries_path = parent_folder + 'geometries/%s' % case

flow = Flow()
flow.load_flow(flow_path, elements_path, geometries_path)

# Prepare pathlines
births   = range(0,1)
lifetime = 50

# Initial values
x0 = -0.7
#list_y = np.array([0.1])
list_y = np.array([0.01, 0.05, 0.1, -0.01, -0.05, -0.1])

initial_positions = prepare_initial_positions(x0, list_y)

#u0 = np.array([10,-5.])
u0 = None

# Arguments to create a list of particles
fill_args = births, lifetime, initial_positions

# Ready to compute
factor = 1

pathlines = spread_pathlines(fill_args, flow, factor, printIt=True)

# Prepare animation
anim_params = {'fps'     :factor,
               'name'    :'../animations/imposedX_atT1/%s_towards_you' % case,
               'fig_save':False}

animation = AnimationPathlines(pathlines, flow, factor, **anim_params)

#animation.elements_on = True
#animation.nodes_on = True
#animation.frames_start = 25*factor
#animation.frames_end = 101*factor
animation.play()

#animation.take_picture(140*factor)
