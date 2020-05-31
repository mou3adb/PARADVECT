"""
This module contains the definition of all linear interpolation functions for a
triangle element.

Their derivatives will serve for the calculation of a variable.

node1, node2, node3 are the coordinates of the three summits of the triangle.

"""
import numpy as np
#==============================================================================
tol = 1e-10

def psi1(*args):
    position, node1, node2, node3 = args

    num   = np.cross(node3-node2, position-node2)
    denom = np.cross(node3-node2, node1-node2)
    rapp  = num/denom

    if abs(rapp) < tol:
        return 0
    else:
        return rapp

def psi1_x(*args):
    position, node1, node2, node3 = args

    x2, y2 = node2
    x3, y3 = node3

    num = y2 - y3
    denom = np.cross(node3-node2, node1-node2)
    return num/denom

def psi1_y(*args):
    position, node1, node2, node3 = args

    x2, y2 = node2
    x3, y3 = node3

    num = x3 - x2
    denom = np.cross(node3-node2, node1-node2)
    return num/denom

def psi2(*args):
    position, node1, node2, node3 = args

    num = np.cross(node1-node3, position-node3)
    denom = np.cross(node3-node2, node1-node2)
    rapp = num/denom

    if abs(rapp) < tol:
        return 0
    else:
        return rapp

def psi2_x(*args):
    position, node1, node2, node3 = args

    x1, y1 = node1
    x3, y3 = node3

    num = y3 - y1
    denom = np.cross(node3-node2, node1-node2)
    return num/denom

def psi2_y(*args):
    position, node1, node2, node3 = args

    x1, y1 = node1
    x3, y3 = node3

    num = x1 - x3
    denom = np.cross(node3-node2, node1-node2)
    return num/denom

def psi3(*args):
    position, node1, node2, node3 = args

    num = np.cross(node2-node1, position-node1)
    denom = np.cross(node3-node2, node1-node2)
    rapp = num/denom

    if abs(rapp) < tol:
        return 0
    else:
        return rapp

def psi3_x(*args):
    position, node1, node2, node3 = args

    x1, y1 = node1
    x2, y2 = node2

    num = y1 - y2
    denom = np.cross(node3-node2, node1-node2)
    return num/denom

def psi3_y(*args):
    position, node1, node2, node3 = args

    x1, y1 = node1
    x2, y2 = node2

    num = x2 - x1
    denom = np.cross(node3-node2, node1-node2)
    return num/denom
#==============================================================================
def interpolation_vector(*args):
    v1 = psi1(*args)
    v2 = psi2(*args)
    v3 = psi3(*args)
    return np.array([v1, v2, v3])

def grad_interpolation_vector(*args):
    v1_x = psi1_x(*args)
    v1_y = psi1_y(*args)
    v2_x = psi2_x(*args)
    v2_y = psi2_y(*args)
    v3_x = psi3_x(*args)
    v3_y = psi3_y(*args)
    return np.array([[v1_x, v1_y],
                     [v2_x, v2_y],
                     [v3_x, v3_y]])