"""
Gmsh format 2.2

"""
import numpy as np

from flow import Flow
from element import Element

from element_search import find_neighbors

from text.text_flow       import write_flow
from text.text_elements   import write_elements
from text.text_geometries import write_geometries
#==============================================================================
def intIt(l):
    return np.array([int(e) for e in l])

def floatIt(l):
    return np.array([float(e) for e in l])

def extract_msh(path_msh):
    f = open(path_msh, 'r')

    nodes_X, nodes_Y = [], []
    elements = []

    line = f.readline()

    # ...
    # $Nodes\n
    # n_nodes
    # ...
    while line != '$Nodes\n':
        line = f.readline()

    line = f.readline()
    n_nodes = int(line.strip())

    for i in range(n_nodes):
        # line = id x y z
        line = f.readline()
        coord = floatIt(line.strip().split())
        nodes_X.append(coord[1])
        nodes_Y.append(coord[2])

    # ...
    # $Elements\n
    # n_elements
    # ...
    while line != '$Elements\n':
        line = f.readline()

    line = f.readline()
    n_elements = int(line.strip())

    count = 0
    for i in range(n_elements):
        # element_id element_type ... ... nodes_id
        line = f.readline()

        coord = intIt(line.strip().split())

        element_type = coord[1]

        if element_type == 9: # 6-node second order triangle
            count += 1
            e = Element(count)
            e.nodes = np.array(coord[-6:])

            elements.append(e)

#        if element_type == 1: # 2-node line
#            e.element_type = 1
#            e.nodes = coord[-2:]
#
#        elif element_type == 2: # 3-node triangle
#            e.element_type = 2
#            e.nodes = coord[-3:]
#
#        elif element_type == 3: # 4-node quadrangle
#            e.element_type = 3
#            e.nodes = coord[-4:]
#
#        elif element_type == 8: # 3-node second order line
#            e.element_type = 8
#            e.nodes = coord[-3:]
#
#        elif element_type == 9: # 6-node second order triangle
#            e.element_type = 9
#            e.nodes = coord[-6:]
#
#        elif element_type == 10: # 9-node second order quadrangle
#            e.element_type = 10
#            e.nodes = coord[-9:]
#
#        elif element_type == 15: # 1-node point
#            e.element_type = 15
#            e.nodes = coord[-1:]
#
#        elements.append(e)

    f.close()

    return np.array(nodes_X), np.array(nodes_Y), np.array(elements)

def generate_poiseuille(path_msh, parent_folder):
    single_nodes_X, single_nodes_Y, elements = extract_msh(path_msh)

    d = np.max(single_nodes_Y) - np.min(single_nodes_Y)
    y_middle = np.min(single_nodes_Y) + d/2

    n_nodes = len(single_nodes_X)

    mu = 1e-3
    p = 2*mu*single_nodes_X

    U = d**2/4 - (single_nodes_Y - y_middle)**2
    V = np.zeros(n_nodes)

    nodes_X, nodes_Y = np.array([]), np.array([])
    Us, Vs, ps = np.array([]), np.array([]), np.array([])

    Nt = 101
    times = np.linspace(0, 1, Nt)
    for t in times:
        nodes_X = np.vstack([nodes_X, single_nodes_X]) if nodes_X.size else single_nodes_X
        nodes_Y = np.vstack([nodes_Y, single_nodes_Y]) if nodes_Y.size else single_nodes_Y

        Us = np.vstack([Us, U]) if Us.size else U
        Vs = np.vstack([Vs, V]) if Vs.size else V
        ps = np.vstack([ps, p]) if ps.size else p

    Re, Ur = 1e-3*1*d/mu, np.inf # Reynolds number and reduced velocity are not
                            # defined in the Hagen-Poiseuille problem

    flow = Flow()

    flow.Re, flow.Ur = Re, Ur

    flow.times = times
    flow.nodes_X, flow.nodes_Y = nodes_X, nodes_Y
    flow.Us, flow.Vs, flow.ps = Us, Vs, ps

    write_flow(flow, parent_folder + 'flows/poiseuille')

    find_neighbors(elements)
    write_elements(elements, parent_folder + 'elements/poiseuille')

    write_geometries(np.array([]), parent_folder + 'geometries/poiseuille')

def generate_periodic(path_msh, parent_folder):
    single_nodes_X, single_nodes_Y, elements = extract_msh(path_msh)

    d = np.max(single_nodes_Y) - np.min(single_nodes_Y)

    Nt = 101
    times = np.linspace(0, 1, Nt)

    period = 0.25
    w = 2*np.pi/period

    # U = U0*cos(wt) with U0 = 1
    # Navier-Stokes, uniform:
    # rho dU/dt + 0 = - dp/dx with rho = 1
    # dp/dx = rhoU0*w*sin(wt)
    # p = p0 + rhoU0*w*sin(wt) with p0 = 0

    nodes_X, nodes_Y = np.array([]), np.array([])
    Us, Vs, ps = np.array([]), np.array([]), np.array([])

    for t in times:
        nodes_X = np.vstack([nodes_X, single_nodes_X]) if nodes_X.size else single_nodes_X
        nodes_Y = np.vstack([nodes_Y, single_nodes_Y]) if nodes_Y.size else single_nodes_Y

        U = 0*nodes_X + np.cos(w*t)
        V = 0*nodes_X
        p = 0*nodes_X + w*np.sin(w*t)

        Us = np.vstack([Us, U]) if Us.size else U
        Vs = np.vstack([Vs, V]) if Vs.size else V
        ps = np.vstack([ps, p]) if ps.size else p

    Re, Ur = 1*1*d/1e-6, np.inf

    flow = Flow()

    flow.Re, flow.Ur = Re, Ur

    flow.times = times
    flow.nodes_X, flow.nodes_Y = nodes_X, nodes_Y
    flow.Us, flow.Vs, flow.ps = Us, Vs, ps

    write_flow(flow, parent_folder + 'flows/periodic')

    find_neighbors(elements)
    write_elements(elements, parent_folder + 'elements/periodic')

    write_geometries(np.array([]), parent_folder + 'geometries/periodic')

def generate_inviscid(path_msh, parent_folder):
    single_nodes_X, single_nodes_Y, elements = extract_msh(path_msh)

    rs     = np.sqrt(single_nodes_X**2 + single_nodes_Y**2)
    thetas = np.arctan2(single_nodes_Y, single_nodes_X)
    
    Ur, Utheta, p = [], [], []
    for r, theta in zip(rs, thetas):
        if r == 0:
            Ur.append(0)
            Utheta.append(0)
            p.append(0)
        else:
            Ur.append((1 - (0.5/r)**2)*np.cos(theta))
            Utheta.append((1 + (0.5/r)**2)*np.sin(theta))
            p.append(2*(0.5/r)**2 * np.cos(2*theta) - (0.5/r)**4)
            
    Ur = np.array(Ur)
    Utheta = np.array(Utheta)
    p = np.array(p)
    
    U = Ur*np.cos(thetas) - Utheta*np.sin(thetas)
    V = Ur*np.sin(thetas) - Utheta*np.cos(thetas)
    
    nodes_X, nodes_Y = np.array([]), np.array([])
    Us, Vs, ps = np.array([]), np.array([]), np.array([])

    Nt = 101
    times = np.linspace(0, 1, Nt)
    for t in times:
        nodes_X = np.vstack([nodes_X, single_nodes_X]) if nodes_X.size else single_nodes_X
        nodes_Y = np.vstack([nodes_Y, single_nodes_Y]) if nodes_Y.size else single_nodes_Y

        Us = np.vstack([Us, U]) if Us.size else U
        Vs = np.vstack([Vs, V]) if Vs.size else V
        ps = np.vstack([ps, p]) if ps.size else p

    Re, Ur = 1e+6, 0.

    flow = Flow()

    flow.Re, flow.Ur = Re, Ur

    flow.times = times
    flow.nodes_X, flow.nodes_Y = nodes_X, nodes_Y
    flow.Us, flow.Vs, flow.ps = Us, Vs, ps

    write_flow(flow, parent_folder + 'flows/potential')

    find_neighbors(elements)
    write_elements(elements, parent_folder + 'elements/potential')

    write_geometries(np.array([[5,407,404,408,405,409,406,410,6,414,411,415,412,416,413,417]]),
                     parent_folder + 'geometries/potential')

