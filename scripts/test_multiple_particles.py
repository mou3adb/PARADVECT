"""
This is a script to test a particle simulation.

Please change any ..._path to your corresponding file path.

"""
import sys
sys.path.append('..')

import numpy as np

from flow import Flow

from animationparticles import AnimationParticles

from text.text_particles import read_particles, write_particles

from test_multiple_particles_scripts import prepare_initial_positions, spread_particles
#==============================================================================
# Prepare flow
parent_folder = r'../'

case = 'uniform'

flow_path       = parent_folder + 'flows/%s' % case
elements_path   = parent_folder + 'elements/%s' % case
geometries_path = parent_folder + 'geometries/%s' % case

#flow = Flow()
#flow.load_flow(flow_path, elements_path, geometries_path)

# Prepare particle
diameter = 0.1
density  = 2.0

births = range(0,1)
lifetime = 20

# Initial values
x0 = -0.75
list_y = np.array([0.])
#list_y = np.array([0.2, 0.423, 0.6])
#list_y = np.array([0.415, 0.6])

initial_positions = prepare_initial_positions(x0, list_y)

#u0 = None
u0 = np.array([0,0.])

# Arguments to create a list of particles
fill_args = diameter, density, births, lifetime, initial_positions, u0

# Ready to compute
factor = 2

particles, captured_ones = spread_particles(fill_args, flow, factor,
                                            printIt=False, too_far_stop=False)

p = particles[0]

# Prepare animation
anim_params = {'fps'     :factor,
               'name'    :'../animations/imposedX_atT1/%s_towards_you' % case,
               'fig_save':False}

#animation = AnimationParticles(particles, flow, factor, **anim_params)
#animation = AnimationParticles(captured_ones, flow, factor, **anim_params)

#animation.paths_on = True
#animation.elements_on = True
#animation.nodes_on = True
#animation.punctual = False
#animation.frames_start = 42*factor
#animation.frames_end = 30*factor
#animation.play()

#animation.take_picture(20*factor)
