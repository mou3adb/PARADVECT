"""
This is a script to test a particle simulation.

Please change any ..._path to your corresponding file path.

"""
import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as pp

from flow import Flow

from animationparticles import AnimationParticles

from text.text_particles import read_particles, write_particles

from test_multiple_particles_scripts import prepare_initial_positions, spread_particles

from axes_world import one_by_two
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
list_factors = np.array([1,2,4,8,16,32,64,128])
#list_trajectories = []
#for factor in list_factors:
#    particles, captured_ones = spread_particles(fill_args, flow, factor,
#                                                printIt=False, too_far_stop=False)
#    
#    traj = particles[0].trajectory
#    list_trajectories.append(traj)
#    
#delta_t = 0.1
#def one_each_two(l):
#    return np.array([l[i] for i in range(0,len(l),2)])
#
#list_epsilon = []
#for i in range(len(list_trajectories)-1):
#    traj_coarse, traj_fine = list_trajectories[i], list_trajectories[i+1]
#    
#    epsilon = np.sqrt( delta_t/list_factors[i] \
#                      * np.sum( (traj_coarse - one_each_two(traj_fine))**2 ) )
#    list_epsilon.append(epsilon)
#list_epsilon_big = np.array(list_epsilon)

# =============================================================================
def plot_discretisation_error_L2(ax, list_dt, errorsL2, color, marker):
    ax.plot(list_dt[1:], errorsL2, linestyle='-',
            linewidth=1, color=color, marker=marker, markeredgecolor='black',
            markeredgewidth=0.5, alpha=0.75)

    ax.set(xscale='log', yscale='log')
    
    ylabel = r'$\varepsilon_{\Delta t}$'
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(r'$\Delta t$', fontsize=12)
    
def plot_estimated_order_L2(ax, list_dt, errorsL2, color, marker):
    obs_p = np.log(errorsL2[:-1]/errorsL2[1:])/np.log(2)
    print(obs_p)
    
    ax.plot(list_dt[2:], obs_p, linestyle='', linewidth=1,
            color=color, marker=marker, markeredgecolor='black',
            markeredgewidth=0.5, alpha=0.75)
    
    ax.set_xscale('log')
    
    ax.set_ylim([0.89, 1.11])
#    ax.set_ylim([1.89, 2.11])
    
    ax.set_yticks([0.9, 0.95, 1, 1.05, 1.1])

    ax.set_xlabel(r'$\Delta t$', fontsize=12)
    ax.set_ylabel(r'$\hat{p}_{\Delta t}$', fontsize=12)
# =============================================================================
list_dt = 0.1/list_factors
ax_a, ax_b = one_by_two()

ax_a.plot(list_dt[1:5], 2e-4*list_dt[:4], linestyle='--',
          color='black')
ax_b.plot(list_dt[2:], 1 + 0*list_dt[2:], linestyle='--',
          color='black')

plot_discretisation_error_L2(ax_a, list_dt, list_epsilon_small, 'blue', '.')
plot_estimated_order_L2(ax_b, list_dt, list_epsilon_small, 'blue', '.')
plot_discretisation_error_L2(ax_a, list_dt, list_epsilon_big, 'gray', 'o')
plot_estimated_order_L2(ax_b, list_dt, list_epsilon_big, 'gray', 'o')
