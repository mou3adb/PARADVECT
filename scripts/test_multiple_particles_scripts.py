"""
These are some scripts used in testing multiple particles.
"""
import time
import datetime

import numpy as np

from particle import Particle
#==============================================================================
def fill_particles(diameter, density, births, lifetime, initial_positions, u0):
    # someParticles = [p, p, ..., p Ninlets times for birth1
    #             then p, p, ..., p Ninlets times for birth2
    #                  ...
    #                  repeated len(births) time]

    someParticles = []
    for b in births:
        list_p = [Particle(diameter, density, b, lifetime, pos0, u0) \
                  for pos0 in initial_positions]

        someParticles.extend(list_p)

    return np.array(someParticles)

def prepare_initial_positions(x0, list_y):
    return np.array([[x0, y] for y in list_y])

def compute_particles(particles, flow, factor, printIt, too_far_stop):
    t1 = time.time()
    cpu_time = 0

    count = 1
    Npar = len(particles)

    for p in particles:
        print('\n---> Particle %d (remaining: %d)' % (count, Npar - count))

        p.compute_trajectory(flow, factor, printIt, too_far_stop)

        cpu_time = time.time() - t1
        print('Accumulated CPU_TIME = %.2f seconds = %s (hh:mm:ss)' \
              % (cpu_time, datetime.timedelta(seconds=cpu_time)))

        count += 1

    captured_ones = np.array([p for p in particles if p.captured])

    return particles, captured_ones

def spread_particles(fill_args, flow, factor, printIt, too_far_stop):
    # fill_args = diameter, density, births, lifetime, initial_positions, u0
    particles = fill_particles(*fill_args)

    return compute_particles(particles, flow, factor, printIt, too_far_stop)