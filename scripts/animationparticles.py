import numpy as np
import matplotlib.pyplot as pp

from matplotlib import rc, rcParams
from matplotlib.animation import FuncAnimation
#==============================================================================
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':'Computer Modern Roman'})
rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

class AnimationParticles(FuncAnimation):
    """
    AnimationParticles class for animation of particles in a flow.

    Attributes
    ----------

    punctual: boolean
        If False, shows the particle as a small circle of dimater equal to
        particle.diameter.

    nodes_on: boolean
        Toogle on/off the nodes of the grid.

    element_on: boolean
        Highlights the element where the particle is laying.

    paths_on: boolean
        Shows the particle trajectory line.

    fig_save: boolean
        Save the image sequence of the animation.

    format: string
        Image format.

    name: string
        Image sequence name.

    fps: integer
        Frames per second, that is, the frequency of image saving
        e.g. fps = 5 --> save one image each 5 frames.

    factor: integer
        Refinement factor of the flow timeline.

    particles: list of Particle objects.

    flow: Flow object

    """

    def __init__(self, particles, flow, factor, fps, name, fig_save):
        self.punctual = True
        self.nodes_on = False
        self.elements_on = False
        self.paths_on = False

        self.fig_save = fig_save
        self.format = 'png'

        self.name = name

        self.fps = fps

        self.factor = factor

        self.particles = particles

        self.flow = flow

        diameter = particles[0].diameter
        self.end_title = r'\text{potential flow}, R = %s' % diameter
#        self.end_title = r'%s, R = %s' % (flow, diameter)

        self.times     = flow.times
        self.new_times = np.linspace(flow.times[0], flow.times[-1],
                                     1 + (len(flow.times)-1)*factor)

        self.dt = flow.times[1] - flow.times[0]

        self.nodes_X = flow.nodes_X
        self.nodes_Y = flow.nodes_Y

        self.geometries = flow.geometries

        # Define the figure and axis
        self.figure = pp.figure(figsize=(5,3.5))
        self.ax = self.figure.add_axes([0.13,0.11,0.85,0.83])
        self.ax.set_xlim([-2,6])
        self.ax.set_ylim([-2,2])
        self.ax.set_aspect(aspect='equal')

        self.xmin, self.xmax = np.min(flow.nodes_X[0]), np.max(flow.nodes_X[0])
        self.ymin, self.ymax = np.min(flow.nodes_Y[0]), np.max(flow.nodes_Y[0])

        self.ax.set_xlabel(r'$X/D$', fontsize=12)
        self.ax.set_ylabel(r'$Y/D$', fontsize=12)

        # Define lines
        self.line_nodes, = self.ax.plot([], [],
                                        marker='.',
                                        color='black',
                                        linestyle='')

        self.lines_elems = self.create_lines_for_elements()

        self.lines_geom = self.create_lines_for_geometries()

        # Frame range of the animation
        self.frames_start = 0
        self.frames_end   = len(self.new_times)

        # To plot circular particles
        self.h = np.linspace(0, 1, 20)

    def create_lines_for_elements(self):
        lines_elements = []

        for i in range(len(self.particles)):
            line, = self.ax.plot([], [],
                                 linestyle='-',
                                 linewidth=1,
                                 color='magenta')
            lines_elements.append(line)

        return lines_elements

    def create_lines_for_geometries(self):
        lines_geometries = []

        for i in range(len(self.geometries)):
            line, = self.ax.plot([], [],
                                 color='black',
                                 linestyle='-',
                                 linewidth=1)
            lines_geometries.append(line)

        return lines_geometries

    def create_lines_for_particles(self):
        lines_particles  = []
        continuous_lines = []

        plot_params = {'color'    :None,
                       'linestyle':None,
                       'marker'   :None,
                       'linewidth':1}

        for p in self.particles:
            if self.punctual:
                plot_params['linestyle'] = ''
                plot_params['marker'] = '.'

                # for fancy pictures
#                plot_params['marker'] = 'o'
#                plot_params['markersize'] = 13
#                plot_params['markeredgecolor'] = 'black'

                if p.captured:
                    plot_params['color'] = 'blue'

                    # for fancy pictures
#                    plot_params['color'] = 'seagreen'
                else:
                    plot_params['color'] = 'red'

                    # for fancy pictures
#                    plot_params['color'] = 'seagreen'

            else:
                plot_params['linestyle'] = '-'
                plot_params['marker'] = ''

                if p.captured:
                    plot_params['color'] = 'blue'
                else:
                    plot_params['color'] = 'red'

            line, = self.ax.plot([], [], **plot_params)

            lines_particles.append(line)

            plot_params['linestyle'] = '-'
            plot_params['marker'] = ''

            continuous_line, = self.ax.plot([], [], **plot_params)

            continuous_lines.append(continuous_line)

        return lines_particles, continuous_lines

    def play(self):
        """
        This function starts the animation.
        """
        # We initiate lines and continuous line here, to avoid the default
        # attribute punctual = True.
        self.lines, self.continuous_lines = self.create_lines_for_particles()

        maximize = False
        if maximize:
            manager = pp.get_current_fig_manager()
            manager.window.showMaximized()

        frames = np.arange(self.frames_start, self.frames_end, self.fps)

        FuncAnimation.__init__(self, self.figure, self.animate,
                               frames=frames, interval=100)

    def take_picture(self, k):
        # Same comment as in self.play()
        self.lines, self.continuous_lines = self.create_lines_for_particles()

        self.animate(k)

        pp.show()

    def remove_ornaments_and_adjust(self, xmin, xmax, ymin, ymax):
        self.ax.axis('off')
        self.ax.set_title('')

        self.ax.set_xlim([xmin,xmax])
        self.ax.set_ylim([ymin,ymax])

    def animate(self, k):
        tk = self.new_times[k]

        pp.title(r'$t = %.2f, %s$' % (tk, self.end_title), fontsize=12)

        # Useful when taking pictures
#        self.ax.plot([], [],
#                     linestyle='',
#                     label=r'$t = %.2f$' % (tk-400))
#        self.ax.legend(loc='upper left',
#                       fontsize=12,
#                       frameon=False,
#                       ncol=1,
#                       labelspacing=0.3,
#                       handlelength=0.)


        # If we want to draw the geometry, elements, or simply visualize nodes,
        # we need to interpolate (in time) their positions.
        n = k//self.factor
        n = n if n < len(self.times)-1 else n-1

        tn   = self.times[n]
        tnp1 = self.times[n+1]

        interp_nodes_X = self.nodes_X[n+1]*(tk - tn)/self.dt \
                       + self.nodes_X[n]*(tnp1 - tk)/self.dt

        interp_nodes_Y = self.nodes_Y[n+1]*(tk - tn)/self.dt \
                       + self.nodes_Y[n]*(tnp1 - tk)/self.dt

        allNodes = np.array(list(zip(interp_nodes_X, interp_nodes_Y)))

        if self.nodes_on:
            self.line_nodes.set_data(interp_nodes_X, interp_nodes_Y)

        # We draw first the particles (and their corresponding elements if on).
        for i in range(len(self.particles)):
            R = self.particles[i].diameter

            particle = self.particles[i]

            # k in FuncAnimation takes the values of frames, which starts from
            # self.frames_start.
            #
            # If birth = 0, the particle is lauched as the animation start.
            # If birth > 0, it waits in its initial position, then is launched.
            real_k = k - particle.birth*self.factor

            if 0 <= real_k/self.factor*1. <= particle.lifetime - 1:
                x_data = particle.trajectory[real_k, 0]
                y_data = particle.trajectory[real_k, 1]

                continuous_x_data = particle.trajectory[:real_k+1, 0]
                continuous_y_data = particle.trajectory[:real_k+1, 1]

            elif real_k/self.factor*1. < 0:
                x_data = particle.trajectory[0, 0]
                y_data = particle.trajectory[0, 1]

                continuous_x_data = particle.trajectory[0:1, 0]
                continuous_y_data = particle.trajectory[0:1, 1]

            else:
                x_data = []
                y_data = []

                continuous_x_data = []
                continuous_y_data = []

            if not(self.punctual) and x_data:
                x_data = x_data + 0.5*R*np.cos(2*np.pi*self.h)
                y_data = y_data + 0.5*R*np.sin(2*np.pi*self.h)

            self.lines[i].set_data(x_data, y_data)

            if self.paths_on:
                self.continuous_lines[i].set_data(continuous_x_data,
                                                  continuous_y_data)

            if self.elements_on:
                elements = particle.elements
                if real_k < 0:
                    e = elements[0]
                    nodes_coords = allNodes[e.nodes - 1]

                    elem_x = nodes_coords[[0,1,2,0],0]
                    elem_y = nodes_coords[[0,1,2,0],1]

                elif 0 <= real_k <= len(elements) - 1:
                    e = elements[real_k]
                    nodes_coords = allNodes[e.nodes - 1]

                    elem_x = nodes_coords[[0,1,2,0],0]
                    elem_y = nodes_coords[[0,1,2,0],1]

                else:
                    elem_x = []
                    elem_y = []

                self.lines_elems[i].set_data(elem_x,elem_y)

        # Draw geometries
        for i in range(len(self.geometries)):
            # Recall that geometry = [nodes that form it]
            geometry = self.geometries[i]

            nodes_coords = allNodes[geometry - 1]
            loop = np.concatenate([nodes_coords, [nodes_coords[0]]])

            self.lines_geom[i].set_data(loop[:,0], loop[:,1])

        if self.fig_save:
            self.figure.savefig(self.name + '_%d.%s' % (k, self.format))
