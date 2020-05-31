"""
Arguments that are in this module
---------------------------------

geometry: list of integers.
    Contains the ids of the nodes forming the geometry.
    e.g. geometry = [node15, node623, ...]

geometries: list of 'geometry's

allNodes: N x 2 numpy array of integers, where N is the number of node in a grid.
    Contains the coordinates of nodes in the domain
    e.g. allNodes = [[node0_x, node0_y], [node1_x, node1_y], ...]
    and  allNodes[0] = [[node0_x, node0_y]]

R: float
    Diameter --- not radius --- of the particle. Take R = 0 if your are dealing
    with punctual particles.

NOTE: node ids start from 1 in this code. This the reason why we substract 1
to the indices, so we can plug them into numpy arrays.

Example of allNodes variable:
allNodes = np.array(list(zip(flow.nodes_X[n], flow.nodes_Y[n]))), where n is your
time step.

"""
import numpy as np
import matplotlib.pyplot as pp
#==============================================================================
def purify_geometry(geometry):
    """
    It's just because CaDyF contains doublons of geometry nodes, when the
    geometry is defined with multiple curves.

    """
    purified = [geometry[0]]

    for i in range(1, len(geometry)):
        node_i = geometry[i]
        if node_i != purified[-1]:
            purified.append(node_i)

    return np.array(purified)

def is_in_geometry(position, R, geometry, allNodes):
    """
    This function checks if a particle of radius R is inside, or even touches,
    the geometry.

    """
    purified = purify_geometry(geometry)

    # Here we substract '-1' to the ids of the geometry nodes because Python
    # list indices start from 0.
    nodes_coord = allNodes[purified[:-1] - 1]
    n_nodes = len(nodes_coord)

    inside_big = []
    in_rounded_corner = None

    for i in range(n_nodes):
        p0, p1, p2 = nodes_coord[[(i-1)%n_nodes, i, (i+1)%n_nodes]]

        p1pos = position - p1
        d1pos = np.linalg.norm(p1pos)

        p1p2 = p2 - p1
        d12  = np.linalg.norm(p1p2)

        p0p1 = p1 - p0
        d01  = np.linalg.norm(p0p1)

        t_p1p2 = p1p2/d12
        n_p1p2 = np.array([t_p1p2[1], -t_p1p2[0]])

        t_p0p1 = p0p1/d01
        n_p0p1 = np.array([t_p0p1[1], -t_p0p1[0]])

        h1 = np.dot(p1pos, n_p1p2)

        sharp_condition = (h1 <= 0.5*R)
        inside_big.append(sharp_condition)

        within_corner = np.cross(n_p0p1, p1pos)>0 and np.cross(p1pos, n_p1p2)>0
        if within_corner:
            if d1pos <= 0.5*R:
                in_rounded_corner = True
            else:
                in_rounded_corner = False

    # If there is at least one False, this means that position lays outside the
    # sharp extended geometry.
    #
    # The second condition accounts for the rounded corners. The particle centre
    # should be inside it (imagine a vertex with a sharp angle).
    if inside_big.count(False) > 0:
        return False

    # Here we are inside the sharp extended geometry. We need to verifiy
    # if the point is in a sharp corner.
    if in_rounded_corner == None:
       return True
    else:
       return in_rounded_corner

def is_in_geometries(position, R, geometries_nodes, allNodes):
    N_geometries = len(geometries_nodes)

    for i in range(N_geometries):
        geometry_nodes = geometries_nodes[i]

        if is_in_geometry(position, R, geometry_nodes, allNodes):
            return True, i

    return False, None

def is_up(position, R, geometry, allNodes):
    """
    If you divide your geometry by two, horizontally, this function verifies
    whether the particle is in the upper or lower side.
    """
    nodes_coord = allNodes[geometry - 1]

    xmin, xmax = np.min(nodes_coord[:,0]), np.max(nodes_coord[:,0])

    barycenter = np.mean(nodes_coord, axis=0)

    bary_pos = position[0] - barycenter[0]

    flag = None

    if xmin - 0.5*R <= position[0] <= xmax + 0.5*R:
        if position[1] >= barycenter[1]:
            flag = True
        else:
            flag = False

    return flag, bary_pos

def get_barycenter(geometry, allNodes):
    purified = purify_geometry(geometry)
    nodes_coord = allNodes[purified - 1]

    return np.mean(nodes_coord, axis=0)

def draw_geometry(geometry, allNodes, ax):
    nodes_coord = allNodes[geometry - 1]

    barycenter = np.mean(nodes_coord, axis=0)

    ax.plot(nodes_coord[:,0], nodes_coord[:,1],
            linestyle='-',
            color='black')

#    ax.plot([barycenter[0]], [barycenter[1]],
#            marker='.',
#            color='black')

    ax.axis('equal')

def draw_geometries(geometries_nodes, allNodes, ax):
    for geometry_nodes in geometries_nodes:
        draw_geometry(geometry_nodes, allNodes, ax)

def verify_inclusion(position, R, geometries_nodes, allNodes, ax):
    """
    Quick function, with drawing, to verify if a particle is inside or not.
    """
    h = np.linspace(0, 1, 100)
    ax.plot(position[0] + 0.5*R*np.cos(2*np.pi*h),
            position[1] + 0.5*R*np.sin(2*np.pi*h),
            color='cyan', linestyle='solid')
    ax.plot(position[0], position[1], marker='.', color='cyan')

    draw_geometries(geometries_nodes, allNodes, ax)

    inside = is_in_geometries(position, R, geometries_nodes, allNodes)
    print('Q: is it inside?\nA: %s' % inside[0])
    print('Q: in which geometry?\nA: %s' % inside[1])

#fig = pp.figure()
#ax = fig.add_axes([0.14,0.13,0.83,0.81])

#allNodes = np.array(zip(flow.nodes_X[0,:], flow.nodes_Y[0,:]))
#geometries_nodes = flow.geometries
#draw_geometries(geometries_nodes, allNodes, ax)
#geometry = flow.geometries[0]
#purified = purify_geometry(geometry)
#draw_geometry(purified, allNodes, ax)
#verify_inclusion(np.array([0.10001,0.5]), 0.1, geometries_nodes, allNodes, ax)

