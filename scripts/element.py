import numpy as np
#==============================================================================
class Element(object):
    """
    This class represents an element in a grid.

    Attributes
    ----------
    id: integer >= 1
        The number identifying the element.

    nodes: numpy array of integers representing the ids of the nodes.

        These are the nodes constituting the element.
        For example, if it's a Tria6, the three first ids are the corners,
        and the rest are the points in the middle of edges.

        NOTE: in my case the CFD code returned me nodes ids >= 1. Please replace
        'nodes - 1' by 'nodes' in all the code if your CFD code returns ids >= 0.

    neighbors: numpy array of integers
        List of ids of the neighboring elements. By definition, a neighbor
        shares an edge with the element.
        An element can have only two neighbors if it is laying on a boundary.

    """
    def __init__(self, Id):
        # All ids start from 1
        self.id = Id

        # Contains node ids
        self.nodes = np.array([], dtype='int32')

        # Neighbors are elements that share an edge with the current one
        # This attribute contains neighboring element ids.
        self.neighbors = np.array([], dtype='int32')
