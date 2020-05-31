"""
This module contains three functions:

F: the applied force on the particle.

compute_Uf: interpolates the value of the fluid speed at the particle position.

compute_gradp: interpolates the value of nabla p at the particle position.
"""
import numpy as np

import interpolation_linear, interpolation_quadratic
#==============================================================================
def F(quadV, Uf, gradp, af, R, rhoplus, Re):
#    r  = quadV[:2]
    up = quadV[2:]

    Rep = Re*R*np.linalg.norm(up - Uf)

    Fd = -0.75*(1./rhoplus)*(1/R)*(up - Uf)*(24./(Re*R))*(1 + 0.15*Rep**0.687)
    Fp = -(1./rhoplus)*gradp

    # Added mass effect with variable, directional mass coefficient
#    theta = np.arctan2(r[1],r[0])

#    u_r     = up[0]*np.cos(theta) + up[1]*np.sin(theta)
#    u_theta =-up[0]*np.sin(theta) + up[1]*np.cos(theta)

#    p_r     = np.abs(u_r)/np.linalg.norm(up)
#    p_theta = np.abs(u_theta)/np.linalg.norm(up)

#    Cm = 0.5 + ((3./8)*p_r + (3./16)*p_theta)* (R/(np.linalg.norm(r)-0.5))**3

    Cm = 0.5

    Fa = (Cm/rhoplus)*af

    dupdt = (Fd + Fp + Fa)/(1 + Cm/rhoplus)

    return np.array([up[0], up[1], dupdt[0], dupdt[1]])

def compute_Uf(position, allUf, allNodes, element):
    points = allNodes[element.nodes - 1][[0, 1, 2]]

    interpol = interpolation_quadratic.interpolation_vector(position, *points)

    pointsUf = allUf[element.nodes - 1]

    return sum(np.transpose(np.transpose(pointsUf)*interpol))

def compute_gradp(position, allps, allNodes, element):
    points = allNodes[element.nodes - 1][[0, 1, 2]]

    gradInterpol = interpolation_linear.grad_interpolation_vector(position, *points)

    pointsp = allps[element.nodes - 1][[0,1,2]]

    return sum(np.transpose(np.transpose(gradInterpol)*pointsp))
