import numpy as np
import matplotlib.pyplot as pp

from scipy.interpolate import splrep, splev

from geometries       import *
from element_search   import find_element, find_element_partrack
from particle_scripts import compute_Uf

#==============================================================================
def boundary_layer(flow, n, theta, dr, id_geometry=0):
    """
    This function calculates the boundary layer thickness for a given angle.

    Attributes
    ----------

    flow: Flow object
        The flow you want to calculate its boundary layer.

    n: integer
        The time step when you want to calculate the boundary layer.

    theta: float
        The anglular position where you want to calculate the thickness.
        Note: theta is in RADIANS!!!

    dr: float
        The radial distance increment starting from the edge of the body.

    """
    elements = flow.elements
    geometry = flow.geometries[id_geometry]
    purified = purify_geometry(geometry)

    allNodes = np.array(list(zip(flow.nodes_X[n], flow.nodes_Y[n])))
    allUf    = np.array(list(zip(flow.Us[n], flow.Vs[n])))

    barycenter = get_barycenter(geometry, allNodes)

    nodes_coord = allNodes[purified[:-1] - 1]

    n_nodes = len(nodes_coord)

    # Looking for the arc which crosses (O, theta)
    t_OM = np.array([-np.cos(theta), np.sin(theta)])
    for i in range(n_nodes):
        p2, p3 = nodes_coord[[i, (i+1)%n_nodes]]

        Op2 = p2 - barycenter
        Op3 = p3 - barycenter

        if np.cross(Op2, t_OM)>0 and np.cross(t_OM, Op3)>0:
            break

    # Intersection point M
    p2p3 = p3 - p2
    d23  = np.linalg.norm(p2p3)

    # We calculate sines to use them in the law of sines
    sin2 = np.cross(Op2, t_OM)/np.linalg.norm(Op2)/np.linalg.norm(t_OM)
    sin3 = np.cross(t_OM, Op3)/np.linalg.norm(Op3)/np.linalg.norm(t_OM)

    # d2 + d3 = d23
    # and
    # sin2/d2 = sin3/d3
    # yield
    # ||
    # ||
    # v
    d2 = d23/(1 + sin2/sin3)

    t_p2p3 = p2p3/d23
    n_p2p3 = np.array([t_p2p3[1], -t_p2p3[0]])

    M = barycenter + Op2 + d2*t_p2p3

    # Pick any faraway point to calculate U_theta (and not U!!)
    # U_theta or -U_theta doesn't change anything, because we take the norm
    faraway_point   = M + 25.*n_p2p3
    faraway_element = find_element(faraway_point, elements, allNodes)

    Ufaraway = compute_Uf(faraway_point, allUf, allNodes, faraway_element)

    norm_Ufaraway_theta = Ufaraway[0]*np.sin(theta) \
                        + Ufaraway[1]*np.cos(theta)

    # Initializing the iterations
    count = 0
    delta = dr
    point = M + delta*n_p2p3

    element = find_element(point, elements, allNodes)

    U = compute_Uf(point, allUf, allNodes, element)

    norm_U_theta = U[0]*np.sin(theta) + U[1]*np.cos(theta)

    while norm_U_theta/norm_Ufaraway_theta < 0.99:
        count += 1

        delta += dr
        point = M + delta*n_p2p3

        element = find_element_partrack(point, element, elements, allNodes)

        U = compute_Uf(point, allUf, allNodes, element)

        norm_U_theta = U[0]*np.sin(theta) + U[1]*np.cos(theta)

    print('Took %d while loops to find delta.' % count)
    print('delta = %e' % delta)

    return delta

#==============================================================================
# Quick script to calculate the thickness for a range of angles
n  = 0
dr = 1e-3

thetas      = np.linspace((5./180), (70./180)*np.pi, 50)
long_thetas = np.linspace((5./180), (70./180)*np.pi, 200)

def calculate_delta(flow, n, thetas, dr):
    deltas = np.empty(0)

    for theta in thetas:
        delta = boundary_layer(flow, n, theta, dr)
        
        deltas = np.append(deltas, delta)

    return deltas

def for_plotting(deltas, thetas, long_thetas):    
    repDeltas    = splrep(thetas, deltas)
    splineDeltas = splev(long_thetas, repDeltas)

    radius = 0.5
    x_deltas = -(radius + splineDeltas)*np.cos(long_thetas)
    y_deltas =  (radius + splineDeltas)*np.sin(long_thetas)
    
    return x_deltas, y_deltas

