"""
This is a script to test a pathline simulation.

Please change any ..._path to your corresponding file path.

"""
import numpy as np

from flow import Flow

from frontline import Frontline

from animationfrontlines import AnimationFrontlines
#==============================================================================
# Prepare flow
parent_folder = r'../'

case = 'potential'

flow_path       = parent_folder + 'flows/%s' % case
elements_path   = parent_folder + 'elements/%s' % case
geometries_path = parent_folder + 'geometries/%s' % case

#flow = Flow()
#flow.load_flow(flow_path, elements_path, geometries_path)

# Prepare pathline
birth    = 0
lifetime = 50

# Initial values
len_segment = 5
x0 = -1.
x0s = np.linspace(  x0,  x0, len_segment)
y0s = np.linspace(-1.5, 1.5, len_segment)
segment = np.array(list(zip(x0s,y0s)))

# Create the frontline
#frontline = Frontline(segment, birth, lifetime)

# Ready to compute
factor = 1

#frontline.compute_frontline(flow, factor, printIt=True)

# Prepare animation
anim_params = {'fps'     :factor,
               'name'    :'../animations/imposedX_atT1/%s_towards_you' % case,
               'fig_save':False}

animation = AnimationFrontlines([frontline], flow, factor, **anim_params)
#animation.nodes_on = True
#animation.frames_start = 41*factor
#animation.frames_end = 161*factor
animation.play()
