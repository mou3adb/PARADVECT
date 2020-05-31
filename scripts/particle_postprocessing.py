"""
This module contains functions used in the post-processing of particle simulation
results.
"""
import numpy as np

import matplotlib.pyplot as pp

#==============================================================================
# Forces calculation
def drag(particle):
    """
    Calculates the drag of the particle. The drag coefficient used here is
    the Schiller-Nauman interpolation.

    """
    R       = particle.diameter
    rhoplus = particle.density
    Re      = particle.Re

    up = particle.velocities
    Uf = particle.fluid_velocities
    relative_u = up - Uf

    Rep = Re*R*np.linalg.norm(relative_u, axis=1)

    # A multiplication with the diagonal of Rep is mathematically correct,
    # but would give an array full of nans.
    column1 = relative_u[:,0] * (1 + 0.15*Rep**0.687)
    column2 = relative_u[:,1] * (1 + 0.15*Rep**0.687)

    return -0.75*(1./rhoplus)*(24./(Re*R**2))*np.array(list(zip(column1, column2)))

def drags(particles):
    return np.array([drag(p) for p in particles])

def pressure(particle):
    """
    Calculates the pressure gradient force, aka Froude-Krylov force.
    """
    rhoplus = particle.density

    gradp = particle.pressure_gradients

    return -(1./rhoplus)*gradp

def pressures(particles):
    return np.array([pressure(p) for p in particles])

def inertia(particle):
    """
    Calculates the added mass force.
    """
    rhoplus = particle.density

    ap = particle.accelerations
    af = particle.fluid_accelerations
    Cm = 0.5

    return -(Cm/rhoplus)*(ap - af)

def inertias(particles):
    return np.array([inertia(p) for p in particles])

#==============================================================================
# Simple extractions
def velocities(particle):
    up = particle.velocities
    Uf = particle.fluid_velocities

    relative = up - Uf

    return up, Uf, relative

def velocity_norms(particle):
    up, Uf, relative = velocities(particle)

    up_norm = np.linalg.norm(up, axis=1)
    Uf_norm = np.linalg.norm(Uf, axis=1)

    relative_norm = np.linalg.norm(relative, axis=1)

    return up_norm, Uf_norm, relative_norm

def accelerations(particle):
    af = particle.fluid_accelerations
    ap = particle.accelerations

    relative = ap - af

    return ap, af, relative

def acceleration_norms(particle):
    ap, af, relative = accelerations(particle)

    ap_norm = np.linalg.norm(ap, axis=1)
    af_norm = np.linalg.norm(af, axis=1)

    relative_norm = np.linalg.norm(relative, axis=1)

    return ap_norm, af_norm, relative_norm

def particle_Reynolds(particle):
    R  = particle.diameter
    Re = particle.Re

    relative = particle.velocities - particle.fluid_velocities

    return Re*R*np.linalg.norm(relative, axis=1)

#def particle_Cm(particle):
#    R = particle.diameter
#
#    r = particle.trajectory
#    up = particle.velocities
#
#    theta = np.arctan2(r[:,1],r[:,0])
#
#    u_r     = up[:,0]*np.cos(theta) + up[:,1]*np.sin(theta)
#    u_theta =-up[:,0]*np.sin(theta) + up[:,1]*np.cos(theta)
#
#    p_r     = np.abs(u_r)/np.linalg.norm(up, axis=1)
#    p_theta = np.abs(u_theta)/np.linalg.norm(up, axis=1)
#
#    Cm = 0.5 + ((3./8)*p_r + (3./16)*p_theta)* (R/(np.linalg.norm(r,axis=1)-0.5))**3
#
#    return Cm

def give_factor(particle):
    # Returns the refinement factor of the timeline.
    # Recall that:
    # If the original timelime is [0, 0.1, ..., 1, 1.1, ...]
    # and that factor = 10, then the new timeline is
    # [0, 0.01, ..., 0.09, 0.1, 0.11, ..., 0.99, 1, 1.01, ...]
    # From this example it is possible to infer the refinement factor formula.

    return len(particle.trajectory)/particle.lifetime

def give_nstop(particle):
    # we assume that finicky = False in particle.py
    up = particle.velocities

    return len(up[np.logical_not(np.isnan(up[:,0]))])

#==============================================================================
# Graphics
def plot_data_vs_time(data, factor, ax, color, label, ls):
#    ax.set_xlabel(r'$\bar{t}$', fontsize=12)
#    ax.set_xlabel(r'$tU_{0}/D$', fontsize=12)

    len_traj = len(data)
    lifetime = len_traj/factor/10.

    t = np.linspace(0,lifetime,len_traj)

#    skip = 10
    skip = 2*factor
    ax.plot(t[::skip], data[::skip], marker='',
            color=color,
            label=label,
            linestyle=ls,
            linewidth=1.5)

#    ax.plot([], [], linestyle='', label=label)
    ax.legend(loc='upper right', numpoints=1,
              fontsize='large',
              markerscale=2,
              frameon=False,
              ncol=1,
              labelspacing=0.5,
              handlelength=1.8)
    
