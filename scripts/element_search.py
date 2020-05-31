import time, datetime

import numpy as np
import matplotlib.pyplot as pp

from interpolation_linear import psi1, psi2, psi3
#==============================================================================
def is_in_element(position, element, allNodes):
    """
    This function tests the inclusion of a point in an element. It calculates
    the form functions. If one if negative, it means that the point is outside.

    NOTE: We substract 1 from element.nodes since node ids start from 1, like
    element ids.

    example of allNodes array:
    allNodes = np.array(list(zip(flow.nodes_X[0],flow.nodes_Y[0])))

    """
    nodes_coords = allNodes[element.nodes - 1]
    p1, p2, p3 = nodes_coords[[0, 1, 2]]

    N1 = psi1(position, p1, p2, p3)
    N2 = psi2(position, p1, p2, p3)
    N3 = psi3(position, p1, p2, p3)

    # It suffices that one interpolation function takes a negative value to
    # deduce that the point is not inside the element.
    if (N1 >= 0) and (N2 >= 0) and (N3 >= 0):
        return True
    else:
        return False

def find_element(position, elements, allNodes):
    for element in elements:
        if is_in_element(position, element, allNodes):
            return element

def find_opposite_neighbor(node, element, elements):
    """
    This function returns the neighbor of 'element' that doesn't contain the
    node 'node'.

    NOTE: Valid only for triangular elements.

    """

    neighbors = elements[element.neighbors - 1]

    for neighbor in neighbors:
        if node not in neighbor.nodes:
            return neighbor

def find_element_partrack(position, current_element, elements, allNodes):
    """
    This function implements the particle tracer algorithm of
    Lohner and Ambrosiano (1990)

    """
    nodes_coords = allNodes[current_element.nodes - 1]
    p1, p2, p3 = nodes_coords[[0, 1, 2]]

    N1 = psi1(position, p1, p2, p3)
    N2 = psi2(position, p1, p2, p3)
    N3 = psi3(position, p1, p2, p3)

    if (N1 >= 0) and (N2 >= 0) and (N3 >= 0):
        # In this case, the particle is still inside 'current_element'.
        return current_element

    else:
        # e.g.
        # If N3 has the most negative value, then index_smallest = 2.
        # Thus point_id_smallest (>= 0) is the farthest node from 'position'.
        index_smallest = np.argmin([N1, N2, N3])
        point_id_smallest = current_element.nodes[index_smallest]

        next_element = find_opposite_neighbor(point_id_smallest,
                                              current_element,
                                              elements)

        return find_element_partrack(position, next_element, elements, allNodes)

def find_neighbors_aux(element, elements):
    """
    This function finds the neighbors of 'element'.

    """
    nodes = element.nodes[[0, 1, 2]]

    for e in elements:
        temp_nodes = e.nodes[[0, 1, 2]]

        # check = [node1_in_or_out?, node2_in_or_out?, node3_in_or_out?]
        check = []
        for node in nodes:
            if node in temp_nodes:
                check.append(True)
            else:
                check.append(False)

        # When an edge is shared by two elements, they have two nodes in
        # common.
        if check.count(True) == 2:
            element.neighbors = np.concatenate([element.neighbors, [e.id]])

def find_neighbors(elements):
    # This function may take a lot of time, since it explores the whole domain,
    # and checks element by element.
    t1 = time.time()
    print('Searching neighbors of every element in the domain.')
    for e in elements:
        find_neighbors_aux(e, elements)

    cpu_time = time.time() - t1
    print('Done!')
    print('CPU_TIME = %.2f SECONDS = %s (HH:MM:SS)' \
          % (cpu_time, datetime.timedelta(seconds=cpu_time)))
#==============================================================================
# Section relating to element drawing
def draw_element(element, allNodes):
    nodes_coords = allNodes[element.nodes - 1]

    print(nodes_coords)

    plot_params = {'linestyle':'-',
                   'color'    :'olive'}

    pp.plot(nodes_coords[[0,1,2,0],0], nodes_coords[[0,1,2,0],1], **plot_params)

def draw_elements(elements, allNodes, color):
    for element in elements:
        draw_element(element, allNodes, color)
