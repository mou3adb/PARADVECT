"""
The outfile structure is the following:

diameter density
birth lifetime
is_captured stuck_to_geometry theta
(blank line)
Re Ur
(blank line)
n_trajectory
x1 y1 up1 vp1 Uf1 Vf1 gradpx1 gradpy1 ap_x1 ap_y1 af_x1 af_y1
x2 y2 up2 vp2 Uf2 Vf2 gradpx2 gradpy2 ap_x2 ap_y2 af_x2 af_y2
...
xNt yNt upNt vpNt UfNt VfNt gradpxNt gradpyNt ap_xN ap_yN af_xN af_yN

"""
import sys
sys.path.append('..')

import numpy as np

from particle import Particle
#==============================================================================
def floatIt(l):
    return np.array([float(e) for e in l])

def intIt(l):
    return np.array([int(e) for e in l])

def write_particle(p, f):
    f.write('%2.3f %1.3f\n' % (p.diameter, p.density))
    f.write('%d %d\n' % (p.birth, p.lifetime))
    f.write('%s %s %s\n' % (p.captured, p.stuck_to_geometry, p.theta))

    f.write('\n') # blank line
    f.write('%d %.1f\n' % (p.Re, p.Ur))

    f.write('\n')
    Nt = len(p.trajectory)
    f.write('%d\n' % Nt)

    for n in range(Nt):
        f.write('%e '*12 % \
                (p.trajectory[n,0],
                 p.trajectory[n,1],
                 p.velocities[n,0],
                 p.velocities[n,1],
                 p.fluid_velocities[n,0],
                 p.fluid_velocities[n,1],
                 p.pressure_gradients[n,0],
                 p.pressure_gradients[n,1],
                 p.accelerations[n,0],
                 p.accelerations[n,1],
                 p.fluid_accelerations[n,0],
                 p.fluid_accelerations[n,1]))

        f.write('\n')

def write_particles(particles, outfile):
    f = open(outfile, 'w')

    Np = len(particles)

    f.write('%d\n' % Np)
    f.write('\n') # blank line

    for p in particles:
        write_particle(p, f)
        f.write('\n')

    f.close()

def read_particle(f, old_version=False):
    # I kept old_version because I had many particles saved before the final
    # update of this function.
    diameter, density = floatIt(f.readline().strip().split())

    birth, lifetime = intIt(f.readline().strip().split())

    if not(old_version):
        str_captured, str_stuck, str_theta = f.readline().strip().split()
        theta = float(str_theta)

    else:
        str_captured, str_stuck = f.readline().strip().split()

    captured = False if str_captured == 'False' else True
    stuck    = None  if str_stuck    == 'None' else int(str_stuck)

    f.readline() # read the blank line
    Re, Ur = floatIt(f.readline().strip().split())

    f.readline()
    Nt = int(f.readline().strip())

    trajectory = []
    velocities = []

    fluid_velocities   = []
    pressure_gradients = []

    accelerations       = []
    fluid_accelerations = []

    for n in range(Nt):
        if old_version:
            x, y, u, v, U, V, gradpx, gradpy \
            = floatIt(f.readline().strip().split())
        else:
            x, y, u, v, U, V, gradpx, gradpy, ap_x, ap_y, af_x, af_y \
            = floatIt(f.readline().strip().split())

        trajectory.append([x, y])
        velocities.append([u, v])

        fluid_velocities.append([U, V])
        pressure_gradients.append([gradpx, gradpy])

        if not(old_version):
            accelerations.append([ap_x, ap_y])
            fluid_accelerations.append([af_x, af_y])

    pos0 = trajectory[0]
    u0 = velocities[0]
    
    p = Particle(diameter, density, birth, lifetime, pos0, u0)

    p.captured, p.stuck_to_geometry = captured, stuck

    p.Re, p.Ur = Re, Ur

    p.trajectory = np.array(trajectory)
    p.velocities = np.array(velocities)

    p.fluid_velocities   = np.array(fluid_velocities)
    p.pressure_gradients = np.array(pressure_gradients)

    if not(old_version):
        p.accelerations       = np.array(accelerations)
        p.fluid_accelerations = np.array(fluid_accelerations)

        p.theta = theta

    return p

def read_particles(infile, old_version=False):
    f = open(infile, 'r')

    Np = int(f.readline())

    f.readline() # read a blank line

    particles = []

    for i in range(Np):
        particles.append(read_particle(f, old_version))

        f.readline()

    f.close()

    return np.array(particles)
