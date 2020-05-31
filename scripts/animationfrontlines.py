import numpy as np
import matplotlib.pyplot as pp

from matplotlib import rc, rcParams
from matplotlib.animation import FuncAnimation
#==============================================================================
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':'Computer Modern Roman'})
rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

class AnimationFrontlines(FuncAnimation):
    def __init__(self, frontlines, flow, factor, fps, name, fig_save):
        self.nodes_on = False

        self.fig_save = fig_save
        self.format = 'png'

        self.name = name

        self.fps = fps

        self.factor = factor

        self.frontlines = frontlines

        self.birth = frontlines[0].birth
        self.lifetime = frontlines[0].lifetime

        self.flow = flow

        self.times = flow.times
        self.new_times = np.linspace(flow.times[0], flow.times[-1],
                                     1 + (len(flow.times)-1)*factor)

        self.dt = flow.times[1] - flow.times[0]

        self.nodes_X = flow.nodes_X
        self.nodes_Y = flow.nodes_Y

        self.geometries = flow.geometries

        # Define the figure and axis
        self.figure = pp.figure()

        self.xmin, self.xmax = np.min(flow.nodes_X[0]), np.max(flow.nodes_X[0])
        self.ymin, self.ymax = np.min(flow.nodes_Y[0]), np.max(flow.nodes_Y[0])
        self.ax = self.figure.add_subplot(111, xlim=(-2,6), ylim=(-2,2))

        self.ax.set_xlabel(r'$X$', size='xx-large')
        self.ax.set_ylabel(r'$Y$', size='xx-large')
        self.figure.gca().set_aspect('equal', adjustable='box')

        # Define lines
        self.lines_fronts = self.create_lines_for_fronts()

        self.line_nodes, = self.ax.plot([], [],
                                        marker='.',
                                        color='black',
                                        linestyle='')

        self.lines_geom = self.create_lines_for_geometries()

        # Frame range of the animation
        self.frames_start = 0
        self.frames_end = len(self.new_times)

    def create_lines_for_geometries(self):
        lines_geometries = []

        for i in range(len(self.geometries)):
            line, = self.ax.plot([], [],
                                 color='black',
                                 linestyle='-')

            lines_geometries.append(line)

        return lines_geometries

    def create_lines_for_fronts(self):
        lines_fronts = []

        for frontline in self.frontlines:
            line, = self.ax.plot([], [],
                                 color='grey',
                                 linestyle='-',
                                 marker='.',
                                 linewidth=1)

            lines_fronts.append(line)

        return lines_fronts

    def fill_gaps(self):
        all_fronts_X, all_fronts_Y = [], []

        new_Nt = len(self.new_times)

        for frontline in self.frontlines:
            true_birth = self.birth*self.factor
            true_death = true_birth + self.lifetime*self.factor + 1

            tmp_X = frontline.front[:,:,0]
            tmp_Y = frontline.front[:,:,1]

            edge_X = np.array([tmp_X[0]]*true_birth)
            edge_Y = np.array([tmp_Y[0]]*true_birth)

            rear_X = np.array([self.xmax]*(new_Nt-true_death))
            rear_Y = np.array([self.ymax]*(new_Nt-true_death))

            mix_trajectory_X = np.concatenate([edge_X, tmp_X, rear_X])
            mix_trajectory_Y = np.concatenate([edge_Y, tmp_Y, rear_Y])

            all_fronts_X.append(mix_trajectory_X)
            all_fronts_Y.append(mix_trajectory_Y)

        return np.array(list(zip(all_fronts_X, all_fronts_Y)))

    def play(self):
        pp.tight_layout()

        frames = np.arange(self.frames_start, self.frames_end, self.fps)

        FuncAnimation.__init__(self, self.figure, self.animate, frames=frames, interval=100)

    def animate(self, k):
        tk = self.new_times[k]

        pp.title(r'$t = %.2f, \text{potential}$' % tk, fontsize='xx-large')

        n = k//self.factor
        n = n if n < len(self.times)-1 else n-1

        tn = self.times[n]
        tnp1 = self.times[n+1]

        interp_nodes_X = self.nodes_X[n+1]*(tk - tn)/self.dt \
                       + self.nodes_X[n]*(tnp1 - tk)/self.dt
        interp_nodes_Y = self.nodes_Y[n+1]*(tk - tn)/self.dt \
                       + self.nodes_Y[n]*(tnp1 - tk)/self.dt

        allNodes = np.array(list(zip(interp_nodes_X, interp_nodes_Y)))

        if self.nodes_on:
            self.line_nodes.set_data(interp_nodes_X, interp_nodes_Y)

        for i in range(len(self.frontlines)):
            frontline = self.frontlines[i]
            real_k = k - frontline.birth*self.factor

            if 0 <= real_k/self.factor*1. <= frontline.lifetime - 1:
                x_data = frontline.front[:, real_k, 0]
                y_data = frontline.front[:, real_k, 1]

            elif real_k/self.factor*1. < 0:
                x_data = frontline.front[:, 0, 0]
                y_data = frontline.front[:, 0, 1]

            else:
                x_data = []
                y_data = []

            self.lines_fronts[i].set_data(x_data, y_data)

        for i in range(len(self.geometries)):
            # Recall that geometry = [nodes that form it]
            geometry = self.geometries[i]

            nodes_coords = allNodes[geometry - 1]
            loop = np.concatenate([nodes_coords, [nodes_coords[0]]])

            self.lines_geom[i].set_data(loop[:,0], loop[:,1])

        if self.fig_save:
            self.figure.savefig(self.name + '_%d.%s' % (k, self.format))
