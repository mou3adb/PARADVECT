from text.text_flow       import read_flow
from text.text_elements   import read_elements
from text.text_geometries import read_geometries
#==============================================================================
class Flow(object):
    """
    Flow is the class that represents the solution of a fluid flow.

    Attributes
    ----------

    Re: integer
        The Reynolds number of the flow.

    Ur: integer, float if necessary
        The reduced velocity, if an object is vibrating.

    times: numpy array of floats.
        The timeline returned by the CFD code.

    nodes_X, nodes_Y: numpy array of floats.
        The x- and y-coordinates of all nodes in the whole domain.

    Us, Vs: numpy array of floats.
        The x- and y-components of the fluid velocity at each node.

    ps: numpy array of floats.
        The pressure at each node.

    elements: numpy array of Element objects.

    geometries: numpy array of 'geometry's, defined by the nodes of their edges.

    """
    def __init__(self):
        attributes = ['Re', 'Ur', 'times',
                      'nodes_X', 'nodes_Y',
                      'Us', 'Vs', 'ps',
                      'elements', 'geometries']

        for attr in attributes:
            setattr(self, attr, None)

    def __str__(self):
        if self.fixed:
            return r'Re = %s, fixed' % self.Re
        else:
            return r'Re = %s, U_{r} = %s' % (self.Re, self.Ur)

    def load_flow(self, flow_path, elements_path, geometries_path):
        """
        This function must be adapted for each CFD code.

        It reads the output of the CFD solution, and assigns it to the
        attributes.

        Your can script your own read_... functions to read the output file of
        your CFD code.
        """
        self.Re, self.Ur, self.times, \
        self.nodes_X, self.nodes_Y, \
        self.Us, self.Vs, self.ps = read_flow(flow_path)

        self.elements = read_elements(elements_path)

        self.geometries = read_geometries(geometries_path)

        # Check if nodes are fixed, in which case the grid is not moving.
        if (self.nodes_X[0] == self.nodes_X[1]).all():
            self.fixed = True

        else:
            self.fixed = False
