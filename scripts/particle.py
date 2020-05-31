"""
"""
import time, datetime
import numpy as np

from particle_scripts import *

from element_search   import find_element, find_element_partrack

from geometries       import is_in_geometries, is_up
#==============================================================================
class Particle(object):
    """
    Particle is the class that represents a particle in a flow.

    Attributes
    ----------

    diameter: float
        Diameter of the particle.

    density: float
        Mass density of the particle.

    birth: integer
        This is relative to the timeline of the flow.
        e.g. If flow.times = [402.0, 402.1, 402.2, ..., 500.0], then birth = 0
        will launch the particle at t = 402.0, and birth = 2 at t = 402.2.

    lifetime: integer
        This is also relative to the timeline.
        e.g. if birth = 0 and lifetime = 2, the life duration of the particle
        is [402.0, 402.1, 402.2].

    pos0, u0: 1 x 2 numpy arrays
        Initial position and velocity of the particle.

        u0 can take the value None, in which case the code will assign it the
        local fluid speed.

    trajectory, velocities, ...: n_traj x 2 numpy array, where n_traj is the number
        of timesteps.

    captured: boolean
        Tells if the particle is captured by any geometry present in the domain.

    struck_to_geometry: integer
        The id of the geometry that the particle hitted.

        NOTE: this version account for only one geometry.

    theta: float
        The angle where the particle is captured.

    is_up: n_traj x 2 numpy array.
        This is a workaround to find out if the particle is above or below the
        geometry. Useful for the function dichotomy in automatic_search_script.

        is_up[i] = [flag, dist].
        flag: boolean
            It's True when the particle is above and False if below.
            If the particle hasn't attained the geometry yet, flag = None.

        dist:float
            It's the x-center-to-center distance between the particle and the
            leftmost point in the geometry.

    Re, Ur = integer, integer (float whenever possible)
        The Reynolds number and reduced velocity of the flow.

    elements: nmupy array of element objects.
        They are the element from which the particles passed through.
        Useful for visuale verification, in animation for instance.

    """
    def __init__(self, diameter, density, birth, lifetime, pos0, u0):
        self.diameter, self.density = diameter, density

        self.birth, self.lifetime = birth, lifetime

        self.pos0 = pos0
        self.u0   = u0

        self.trajectory = np.array([])
        self.velocities = np.array([])

        self.fluid_velocities    = np.array([])
        self.pressure_gradients  = np.array([])

        self.accelerations = np.array([])
        self.fluid_accelerations = np.array([])

        self.captured = False
        self.stuck_to_geometry = None
        self.theta = np.nan

        # This attribute is useful for the automatic search sript, to see if
        # the particle went above or below the object.
        self.is_up = []

        self.Re = None
        self.Ur = None

        # This attribute is useful to verify if the particle was indeed
        # included in elements.
        self.elements = []

    def __str__(self):
        return 'Particle diameter: %s\nInitial position: %s' \
                % (self.diameter, self.pos0)

    def compute_trajectory(self, flow, factor, printIt, too_far_stop=True):
        """
        Arguments
        ---------

        flow: Flow object

        factor: integer
            The refinement factor of the flow timeline.

        printIt: boolean

        too_far_stop: boolean
            is True by default, because it's convenient in automatic search and
            dichotomy. For tests or other simulations, you need to change it to
            False.

        """
        t1 = time.time()

        print('Stop when too far is %s' % too_far_stop)

        if printIt:
            print(self)
            print('Computing particle trajectory (%s)' % datetime.datetime.now())

        # Loading flow variables
        times = flow.times[self.birth:self.birth+self.lifetime+1]

        nodes_X = flow.nodes_X[self.birth:self.birth+self.lifetime+1]
        nodes_Y = flow.nodes_Y[self.birth:self.birth+self.lifetime+1]

        Us = flow.Us[self.birth:self.birth+self.lifetime+1]
        Vs = flow.Vs[self.birth:self.birth+self.lifetime+1]
        ps = flow.ps[self.birth:self.birth+self.lifetime+1]

        elements = flow.elements

        geometries = flow.geometries

        trio = self.diameter, self.density, flow.Re

        self.Re = flow.Re
        self.Ur = flow.Ur

        # This is to calculate the X-barycenter mean position through time, to
        # include it as a loop break condition, in case the particle goes too far.
        #
        # The geomtry tag is put by default to 0, because in this version only
        # one object is accounted for. Furter impovements are required.
        if len(geometries):
            geometry_tag = 0
            geometry_nodes = geometries[geometry_tag]

            # Note that we aren't using axis=0 arugment in np.mean, because we are
            # averaging the center position through time.
            mean_barycenter_X = np.mean(nodes_X[:, geometry_nodes])

        else:
            mean_barycenter_X = np.max(nodes_X) - 1.5

        # Sometimes time subdivision can be irregular. Otherwise dt is constant.
        Nt  = len(times) # Nt == self.lifetime + 1
        dts = np.diff(times)

        # Here we create a refined timeline.
        new_times = np.linspace(times[0], times[-1], 1 + (Nt-1)*factor)

        new_Nt  = len(new_times)
        new_dts = np.diff(new_times)

        # Index 0 means t0, i.e. we haven't started the integration loop yet.
        if printIt:
            print('t = %f' % new_times[0])

        # Assign initial position.
        position = np.copy(self.pos0)
        self.trajectory = np.vstack([np.copy(self.pos0)])

        # Compute pressure and fluid velocity in that position
        # Here zero index means at t0.
        # e.g.
        # allNodes = [[node1_x, node1_y],
        #             [node2_x, node2_y], ...]
        #
        # allUf    = [[U(node1)_x, V(node1)_y],
        #             [U(node2)_x, V(node2)_y], ...]
        #
        # allps    = [p(node1),
        #             p(node2), ...]
        allNodes     = np.array(list(zip(nodes_X[0], nodes_Y[0])))
        allUf, allps = np.array(list(zip(Us[0], Vs[0]))), ps[0]

        # Search for the element containing the particle.
        element = find_element(position, elements, allNodes)
        self.elements.append(element)

        # Sometimes the particle is chosen inside or at the boundary of an
        # obstacle.
        captured, geometry_tag = is_in_geometries(position,
                                                  self.diameter,
                                                  geometries,
                                                  allNodes)
        if captured:
            assert(not(captured)), 'The initial position of the particle is ' \
                                 + 'inside the geometry. Please choose one ' \
                                 + 'outside and try again.'

        # Check if particle is above of below the obstacle.
        #
        # Notice that this version of the code considers only
        # one geometry. Further improvement must be done to include
        # multiple geometries.
        if len(geometries):
            flag, d_from_left_edge = is_up(position,
                                           self.diameter,
                                           geometries[0],
                                           allNodes)

            self.is_up = [[flag, d_from_left_edge]]

        # Interpolate pressure and fluid velocity.
        Uf    = compute_Uf(position, allUf, allNodes, element)
        gradp = compute_gradp(position, allps, allNodes, element)

        self.fluid_velocities   = np.vstack([np.copy(Uf)])
        self.pressure_gradients = np.vstack([np.copy(gradp)])

        # Calculate the local acceleration.
        #
        # Notice that it isn't /new_dts[0] but indeed /dts[0]. The interpolation
        # is linear, so we need to take the slope.
        AccUs = (Us[1] - Us[0])/dts[0]
        AccVs = (Vs[1] - Vs[0])/dts[0]
        allAccf = np.array(list(zip(AccUs, AccVs)))

        Accf = compute_Uf(position, allAccf, allNodes, element)

        self.fluid_accelerations = np.vstack([np.copy(Accf)])

        # Assign initial velocity.
        if self.u0 is None:
            up = np.copy(Uf)

            self.u0 = np.copy(up)
            self.velocities = np.vstack([np.copy(up)])

        else:
            up = np.copy(self.u0)
            self.velocities = np.vstack([np.copy(up)])

        # We will take this zero acceleration out in the end of the loop. This
        # just temporary in order to initialize the self.accelerations array.
        self.accelerations = np.vstack([np.array([np.nan,np.nan])])

        # Begin time integration. Notice that we start from i = 1.
        nstop = new_Nt-1
        for i in range(1, new_Nt):

            # This is just an arbitrary condition. Change it if necessary.
            if too_far_stop and position[0] > 1. + mean_barycenter_X:
                print('The particle is going too far. The integration is stopped.')

                nstop = i
                break

            # We have the variables at i-1. we are ready to integrate the
            # momentum equation to find the actual position and velocity.
            quadV = np.array([position[0], position[1], up[0], up[1]])

            # Basic Forward Euler Scheme.
            dF = new_dts[i-1]*F(quadV, Uf, gradp, Accf, *trio)

            position += dF[:2]
            up       += dF[2:]

            # After computing the new position (and velocity), we need to check
            # if at that point the particle hits the object.

            ti = new_times[i]

            # n is the time index corresponding to the original timeline.
            #
            # e.g.
            # If the original timeline is [0, 0.1, 0.2, ...] and factor is 10,
            # then the new timeline is [0, 0.01, ..., 0.09, 0.1, 0.11, ...].
            # So if i = 5, then n = 0, whereas if i = 10, then n = 1.
            # These will be used for the coming interpolations.
            n = i//factor
            n = n if n < Nt-1 else n-1

            tn      = times[n]
            tnplus1 = times[n+1]
            dt      = dts[n]

            # If not, we can add the new position and velocity to the lists
            # trajectory and velocities.
            self.trajectory = np.vstack([self.trajectory, position])
            self.velocities = np.vstack([self.velocities, up])

            self.accelerations = np.vstack([self.accelerations, dF[2:]/new_dts[i-1]])

            # This is the state of the element between the two instants n and
            # n+1.
            interp_nodes_X = nodes_X[n+1]*(ti - tn)/dt \
                           + nodes_X[n]*(tnplus1 - ti)/dt

            interp_nodes_Y = nodes_Y[n+1]*(ti - tn)/dt \
                           + nodes_Y[n]*(tnplus1 - ti)/dt

            interp_allNodes = np.array(list(zip(interp_nodes_X, interp_nodes_Y)))

            # Check if particle is on top or bottom of the obstacle.
            #
            # Again, this version of the code considers only one geometry.
            if len(geometries):
                flag, d_from_left_edge = is_up(position,
                                               self.diameter,
                                               geometries[0],
                                               interp_allNodes)

                self.is_up.append([flag, d_from_left_edge])

            # Check if particle is inside an obstacle.
            captured, geometry_tag = is_in_geometries(position,
                                                      self.diameter,
                                                      geometries,
                                                      interp_allNodes)

            if captured:
                nstop = i
                self.captured = True
                self.stuck_to_geometry  = geometry_tag

                self.theta = -np.arctan(position[1]/position[0])

                # We don't care about next values, just put np.nan.
                # BTW, this is useful when plotting. The curv ewill end at this
                # point. This is better than saving n_stop as a nattribute and
                # telling it to stop at n_stop.
                self.fluid_velocities     = np.vstack([self.fluid_velocities,
                                                       np.array([np.nan, np.nan])])

                self.pressure_gradients   = np.vstack([self.pressure_gradients,
                                                       np.array([np.nan, np.nan])])

                self.fluid_accelerations  = np.vstack([self.fluid_accelerations,
                                                       np.array([np.nan, np.nan])])
                break

            if printIt:
                print('t = %.6f' % new_times[i])

            # Again, this is the interpolated states between n and n+1.
            interp_Us = Us[n+1]*(ti - tn)/dt + Us[n]*(tnplus1 - ti)/dt
            interp_Vs = Vs[n+1]*(ti - tn)/dt + Vs[n]*(tnplus1 - ti)/dt
            interp_ps = ps[n+1]*(ti - tn)/dt + ps[n]*(tnplus1 - ti)/dt

            interp_accUs = (Us[n+1] - Us[n])/dt
            interp_accVs = (Vs[n+1] - Vs[n])/dt

            # Interpolated acceleration is constant piecewise, because the
            # velocity interpolation is linear.

            # After interpolating, we create arrays of variables.
            allUf, allps = np.array(list(zip(interp_Us, interp_Vs))), interp_ps

            allAccf = np.array(list(zip(interp_accUs, interp_accVs)))

            # Now that we have all ingredients, we can search the new element
            # containing the particle. this is done using the particle tracking
            # algorithm. It needs the previous element element, the actual
            # position and the node coordinates of the domain allNodes.

            element = find_element_partrack(position, element,
                                            elements, interp_allNodes)

            self.elements.append(element)

            # We can now interpolate the values of pressure and fluid velocity.
            Uf    = compute_Uf(position, allUf, interp_allNodes, element)
            gradp = compute_gradp(position, allps, interp_allNodes, element)

            # We can use compute_Uf here. What it does is simply interpolate
            # the values using P2 elements.
            Accf  = compute_Uf(position, allAccf, interp_allNodes, element)

            self.fluid_velocities   = np.vstack([self.fluid_velocities, Uf])
            self.pressure_gradients = np.vstack([self.pressure_gradients, gradp])

            self.fluid_accelerations = np.vstack([self.fluid_accelerations, Accf])

        # We have now left the time integration loop. There have two choices:
        #
        # 1/ either the particle is not captured, then the for loop has been
        # entirely completed, and len(self.trajectory) = new_Nt,
        #
        # 1bis/ we stopped the simulation because the particle went too far
        #
        # 2/ the particle is captured. In this case, we need to fill in the
        # rest of the trajectory and velocities with the position and velocity
        # of the object.
        if self.captured:
            print('>>>>> CAPTURED! <<<<<')

            for i in range(nstop, new_Nt-1):
                if printIt:
                    print('t = %.6f' % new_times[i+1])

                ti     = new_times[i]
                new_dt = new_dts[i]

                # Notice the 'i-1' index as reference
                n = i//int(factor)
                tn      = times[n]
                tnplus1 = times[n+1]
                dt      = dts[n]

                geometry_nodes = geometries[geometry_tag]

                barycenter_Xnp1 = np.mean(nodes_X[n+1, geometry_nodes])
                barycenter_Ynp1 = np.mean(nodes_Y[n+1, geometry_nodes])

                barycenter_Xn = np.mean(nodes_X[n, geometry_nodes])
                barycenter_Yn = np.mean(nodes_Y[n, geometry_nodes])

                last_position = np.copy(self.trajectory[-1])
                last_position += (new_dt/dt) \
                                * np.array([barycenter_Xnp1 - barycenter_Xn,
                                            barycenter_Ynp1 - barycenter_Yn])

                last_velocity = np.copy(self.velocities[-1])
                last_velocity += (1./dt) \
                                * np.array([barycenter_Xnp1 - barycenter_Xn,
                                            barycenter_Ynp1 - barycenter_Yn])

                self.trajectory = np.vstack([self.trajectory, last_position])

                # If you care about the value of the velocity (which is the
                # velocity of the moving object), uncomment this line.
#                self.velocities = np.vstack([self.velocities, last_velocity])

                self.velocities = np.vstack([self.velocities,
                                             np.array([np.nan, np.nan])])

                self.accelerations = np.vstack([self.accelerations,
                                                np.array([np.nan, np.nan])])

                # If you want to obtain the right values of the fluid
                # velocity and pressure gradient at the ensuing positions,
                # just turn finicky to True. However, make sure that the particle,
                # when stuck to the object, has it center within an element of
                # the mesh. Otherwise, you'll get an error.
                #
                # This step can take longer computational time.
                finicky = False

                if finicky:
                    interp_nodes_X = nodes_X[n+1]*(ti - tn)/dt \
                                   + nodes_X[n]*(tnplus1 - ti)/dt

                    interp_nodes_Y = nodes_Y[n+1]*(ti - tn)/dt \
                                   + nodes_Y[n]*(tnplus1 - ti)/dt

                    interp_allNodes = np.array(list(zip(interp_nodes_X, interp_nodes_Y)))

                    allUf, allps = np.array(list(zip(interp_Us, interp_Vs))), interp_ps

                    element = find_element_partrack(position, element,
                                                    elements, interp_allNodes)

                    self.elements.append(element)

                    Uf    = compute_Uf(position, allUf, interp_allNodes, element)
                    gradp = compute_gradp(position, allps, interp_allNodes, element)

                else:
                    Uf    = np.array([np.nan, np.nan])
                    gradp = np.array([np.nan, np.nan])

                    Accf = np.array([np.nan, np.nan])

                self.fluid_velocities   = np.vstack([self.fluid_velocities, Uf])
                self.pressure_gradients = np.vstack([self.pressure_gradients, gradp])

                self.fluid_accelerations = np.vstack([self.fluid_accelerations, Accf])

        # When it isn't captured.
        else:
            max_domain_X, max_domain_Y = np.max(nodes_X), np.max(nodes_Y)

            corner = np.array([max_domain_X, max_domain_Y])

            # This loop goes in only when too_far_stop = True, because it's the
            # case where new_Nt - 1 != 0.
            for i in range(nstop, new_Nt-1):
                self.trajectory = np.vstack([self.trajectory, corner])

                self.velocities = np.vstack([self.velocities,
                                             np.array([np.nan, np.nan])])

                self.accelerations = np.vstack([self.accelerations,
                                                np.array([np.nan, np.nan])])

                self.fluid_velocities    = np.vstack([self.fluid_velocities,
                                                      np.array([np.nan, np.nan])])

                self.pressure_gradients  = np.vstack([self.pressure_gradients,
                                                      np.array([np.nan, np.nan])])

                self.fluid_accelerations = np.vstack([self.fluid_accelerations,
                                                      np.array([np.nan, np.nan])])

        # We remove the first value we used to initialize vstack in numpy.
        self.accelerations = self.accelerations[1:]
        self.accelerations = np.vstack([self.accelerations,
                                        np.array([np.nan, np.nan])])

        self.is_up = np.array(self.is_up, dtype=object)

        d = time.time() - t1
        print('===/ Trajectory successfully computed :)')
        print('===/ CPU_TIME = %.2f seconds = %s (hh:mm:ss)\n' \
              % (d, datetime.timedelta(seconds=d)))
