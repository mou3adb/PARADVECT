"""
This module contains the definition of all quadratic interpolation functions
for a triangle element.

Their derivatives will serve for the calculation of a variable.

*args = position, node1, node2, node3

node1, node2, node3 are the coordinates of the three summits of the triangle.

"""
import numpy as np

from interpolation_linear import psi1  , psi2  , psi3
from interpolation_linear import psi1_x, psi2_x, psi3_x
from interpolation_linear import psi1_y, psi2_y, psi3_y
#==============================================================================
def Psi1(*args):
    L1 = psi1(*args)

    return L1*(2*L1 - 1)

def Psi1_x(*args):
    L1 = psi1(*args)
    L1_x = psi1_x(*args)

    return (4*L1 - 1)*L1_x

def Psi1_y(*args):
    L1 = psi1(*args)
    L1_y = psi1_y(*args)

    return (4*L1 - 1)*L1_y

def Psi2(*args):
    L2 = psi2(*args)

    return L2*(2*L2 - 1)

def Psi2_x(*args):
    L2 = psi2(*args)
    L2_x = psi2_x(*args)

    return (4*L2 - 1)*L2_x

def Psi2_y(*args):
    L2 = psi2(*args)
    L2_y = psi2_y(*args)
    return (4*L2 - 1)*L2_y

def Psi3(*args):
    L3 = psi3(*args)

    return L3*(2*L3 - 1)

def Psi3_x(*args):
    L3 = psi3(*args)
    L3_x = psi3_x(*args)

    return (4*L3 - 1)*L3_x

def Psi3_y(*args):
    L3 = psi3(*args)
    L3_y = psi3_y(*args)

    return (4*L3 - 1)*L3_y

def Psi4(*args):
    L1 = psi1(*args)
    L2 = psi2(*args)

    return 4*L1*L2

def Psi4_x(*args):
    L1 = psi1(*args)
    L2 = psi2(*args)

    L1_x = psi1_x(*args)
    L2_x = psi2_x(*args)

    return 4*(L1_x*L2 + L1*L2_x)

def Psi4_y(*args):
    L1 = psi1(*args)
    L2 = psi2(*args)

    L1_y = psi1_y(*args)
    L2_y = psi2_y(*args)

    return 4*(L1_y*L2 + L1*L2_y)

def Psi5(*args):
    L2 = psi2(*args)
    L3 = psi3(*args)

    return 4*L2*L3

def Psi5_x(*args):
    L2 = psi2(*args)
    L3 = psi3(*args)

    L2_x = psi2_x(*args)
    L3_x = psi3_x(*args)

    return 4*(L2_x*L3 + L2*L3_x)

def Psi5_y(*args):
    L2 = psi2(*args)
    L3 = psi3(*args)

    L2_y = psi2_y(*args)
    L3_y = psi3_y(*args)

    return 4*(L2_y*L3 + L2*L3_y)

def Psi6(*args):
    L1 = psi1(*args)
    L3 = psi3(*args)

    return 4*L3*L1

def Psi6_x(*args):
    L1 = psi1(*args)
    L3 = psi3(*args)

    L1_x = psi1_x(*args)
    L3_x = psi3_x(*args)

    return 4*(L3_x*L1 + L3*L1_x)

def Psi6_y(*args):
    L1 = psi1(*args)
    L3 = psi3(*args)

    L1_y = psi1_y(*args)
    L3_y = psi3_y(*args)

    return 4*(L3_y*L1 + L3*L1_y)
#==============================================================================
def interpolation_vector(*args):
    v1 = Psi1(*args)
    v2 = Psi2(*args)
    v3 = Psi3(*args)
    v4 = Psi4(*args)
    v5 = Psi5(*args)
    v6 = Psi6(*args)

    return np.array([v1, v2, v3, v4, v5, v6])

def grad_interpolation_vector(*args):
    v1_x = Psi1_x(*args)
    v1_y = Psi1_y(*args)
    v2_x = Psi2_x(*args)
    v2_y = Psi2_y(*args)
    v3_x = Psi3_x(*args)
    v3_y = Psi3_y(*args)
    v4_x = Psi4_x(*args)
    v4_y = Psi4_y(*args)
    v5_x = Psi5_x(*args)
    v5_y = Psi5_y(*args)
    v6_x = Psi6_x(*args)
    v6_y = Psi6_y(*args)

    return np.array([[v1_x, v1_y],
                     [v2_x, v2_y],
                     [v3_x, v3_y],
                     [v4_x, v4_y],
                     [v5_x, v5_y],
                     [v6_x, v6_y]])