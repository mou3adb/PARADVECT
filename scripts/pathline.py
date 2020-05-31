import time, datetime
import numpy as np

from particle_scripts import compute_Uf

from element_search import find_element, find_element_partrack
#==============================================================================
class Pathline(object):
    """
    Pathline class calculates the pathline of a flow

    Attributes
    ----------

    pos0: 2 x 1 numpy array
        Initial position of the fluid particle.

    birth: integer
        The time when you start to launch your dye. It is relative to the
        timeline of the flow.
        e.g. birth = 0 is the beginning of the flow.
             birth = int(len(flow.times)/2.) is the middle of the flow

    lifetime: integer
        The lifetime of your dye.
        Should be in the same units of the flow timeline.

    path: n_path x 2 numpy array, where n_path is the number of time steps.
        The trajectory of the fluid particle.

    elements: python list of Element objects (NOT a numpy array)
        The list of the elements from where the fluid particle has passed
        through.

    """
    def __init__(self, birth, lifetime, pos0):
        self.pos0 = pos0

        self.birth = birth
        self.lifetime = lifetime

        self.path = np.array([])

        self.elements = []

    def compute_pathline(self, flow, factor, printIt):
        t1 = time.time()

        # Loading flow variables
        times = flow.times[self.birth:self.birth+self.lifetime+1]

        nodes_X = flow.nodes_X[self.birth:self.birth+self.lifetime+1]
        nodes_Y = flow.nodes_Y[self.birth:self.birth+self.lifetime+1]

        Us = flow.Us[self.birth:self.birth+self.lifetime+1]
        Vs = flow.Vs[self.birth:self.birth+self.lifetime+1]

        elements = flow.elements

        # Sometimes time subdivision can be irregular. Otherwise dt is constant.
        Nt = len(times) # Nt == self.lifetime + 1
        dts = np.diff(times)

        # Here we create a refined timeline.
        new_times = np.linspace(times[0], times[-1], 1 + (len(times)-1)*factor)
        new_N = len(new_times) - 1
        new_dts = np.diff(new_times)

        # Assign initial position.
        position = np.copy(self.pos0)

        self.path = np.array([np.copy(position)])

        # Compute pressure and fluid velocity in that position
        # Here zero index means at t0.
        # e.g.
        # allNodes = [[node1_x, node1_y],
        #             [node2_x, node2_y], ...]
        # allUf    = [[U(node1)_x, V(node1)_y],
        #             [U(node2)_x, V(node2)_y], ...]
        # allps    = [p(node1),
        #             p(node2), ...]
        allNodes = np.array(list(zip(nodes_X[0], nodes_Y[0])))
        allUf = np.array(list(zip(Us[0], Vs[0])))

        # Search for the element containing the particle.
        element = find_element(position, elements, allNodes)
        self.elements.append(element)

        # Interpolate fluid velocity.
        Uf = compute_Uf(position, allUf, allNodes, element)

        # Start time integration. Notice that we start from i = 1.
        for i in range(1, new_N):
            # Basic Forward Euler Scheme.
            position += Uf*new_dts[i-1]

            self.path = np.vstack([self.path, position])

            ti = new_times[i]

            if printIt:
                print('t = %.6f' % new_times[i])

            # n is the time index corresponding to the original timeline.
            # e.g.
            # If the original timline is [0, 0.1, 0.2, ...] and factor is 10,
            # then the new timeline is [0, 0.01, ..., 0.09, 0.1, 0.11, ...].
            # So if i = 5, then n = 0, whereas if i = 10, then n = 1.
            # These will be used for the coming interpolations.
            n = i//factor
            n = n if n < Nt-1 else n-1

            tn = times[n]
            tnplus1 = times[n+1]
            dt = dts[n]

            # This is the state of the element between the two instants n and
            # n+1.
            interp_nodes_X = nodes_X[n+1]*(ti - tn)/dt \
                           + nodes_X[n]*(tnplus1 - ti)/dt
            interp_nodes_Y = nodes_Y[n+1]*(ti - tn)/dt \
                           + nodes_Y[n]*(tnplus1 - ti)/dt

            interp_allNodes = np.array(list(zip(interp_nodes_X, interp_nodes_Y)))

            # Again, this is the interpolated states between n and n+1.
            interp_Us = Us[n+1]*(ti - tn)/dt + Us[n]*(tnplus1 - ti)/dt
            interp_Vs = Vs[n+1]*(ti - tn)/dt + Vs[n]*(tnplus1 - ti)/dt

            # After interpolating, we create arrays of variables.
            allUf = np.array(list(zip(interp_Us, interp_Vs)))

            # Now that we have all ingredients, we can search the new element
            # containing the particle. this is done using the particle tracking
            # algorithm. It needs the previous element element, the actual
            # position and the node coordinates of the domain allNodes.

            element = find_element_partrack(position, element,
                                            elements, interp_allNodes)

            self.elements.append(element)


            # We can now interpolate the values of pressure and fluid velocity.
            Uf = compute_Uf(position, allUf, interp_allNodes, element)

        d = time.time() - t1
        print('===/ Pathline successfully computed :)')
        print('===/ CPU_TIME = %.2f seconds = %s (hh:mm:ss)\n' \
              % (d, datetime.timedelta(seconds=d)))
