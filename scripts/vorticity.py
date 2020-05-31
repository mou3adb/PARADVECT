import numpy as np

from element_search import find_element

from particle_scripts import compute_Uf
#==============================================================================
def vorticity(flow, n, theta_plus, dr, id_geometry=0):
    elements = flow.elements
    geometry_nodes = flow.geometries[id_geometry]

    allNodes = np.array(list(zip(flow.nodes_X[n], flow.nodes_Y[n])))
    allUf = np.array(list(zip(flow.Us[n], flow.Vs[n])))

    nodes_coord = allNodes[geometry_nodes - 1]
    barycenter = np.mean(nodes_coord, axis=0)

    n_nodes = len(geometry_nodes)

    # Looking for the arc which crosses (O, theta_plus)
    for i in range(n_nodes - 1):
        p2, p3 = nodes_coord[[i, i+1]]

        Op2 = p2 - barycenter
        Op3 = p3 - barycenter

        theta_plus2 = np.pi - np.arctan2(Op2[1], Op2[0])
        theta_plus3 = np.pi - np.arctan2(Op3[1], Op3[0])

        if (theta_plus2 - theta_plus) * (theta_plus3 - theta_plus) <= 0:
            break

    # Intersection point M
    p2p3 = p3 - p2
    d23 = np.linalg.norm(p2p3)

    t_p2p3= p2p3/d23
    n_p2p3 = np.array([t_p2p3[1], -t_p2p3[0]])

    theta_plus2 = np.pi - np.arctan2(Op2[1], Op2[0])
    theta_plus3 = np.pi - np.arctan2(Op3[1], Op3[0])

    norm_Op2 = np.linalg.norm(Op2)
    norm_OM = norm_Op2*np.cos(theta_plus2 - theta_plus)

    # Quick verification
    norm_Op3 = np.linalg.norm(Op3)
    norm_OM_via_p3 = norm_Op3*np.cos(theta_plus - theta_plus3)
    if abs(norm_OM_via_p3 - norm_OM) > 1e-3:
        print('norm via p2 = %s' % norm_OM)
        print('norm via p3 = %s' % norm_OM_via_p3)
        raise Exception('The two norms should be equal. Something wrong in the' \
                        + ' algebra.')

    OM = barycenter + norm_OM*np.array([-np.cos(theta_plus), np.sin(theta_plus)])

    nearby_point = OM + dr*n_p2p3
    nearby_element = find_element(nearby_point, elements, allNodes)

    Unearby = compute_Uf(nearby_point, allUf, allNodes, nearby_element)

    norm_Unearby_theta = np.dot(Unearby,
                                np.array([np.sin(theta_plus),
                                          np.cos(theta_plus)]))

    return norm_Unearby_theta/dr
