"""
These are some scripts used in testing multiple pathlines.
"""
import time
import datetime

import numpy as np

from pathline import Pathline
#==============================================================================
def fill_pathlines(births, lifetime, initial_positions):
    # someParticles = [p, p, ..., p Ninlets times for birth1
    #             then p, p, ..., p Ninlets times for birth2
    #                  ...
    #                  repeated len(births) time]

    somePathlines = []
    for b in births:
        list_p = [Pathline(b, lifetime, pos0) \
                  for pos0 in initial_positions]

        somePathlines.extend(list_p)

    return np.array(somePathlines)

def prepare_initial_positions(x0, list_y):
    return np.array([[x0, y] for y in list_y])

def compute_pathlines(pathlines, flow, factor, printIt):
    t1 = time.time()
    cpu_time = 0

    count = 1
    Npar = len(pathlines)

    for p in pathlines:
        print('\n---> Pathline %d (remaining: %d)' % (count, Npar - count))

        p.compute_pathline(flow, factor, printIt)

        cpu_time = time.time() - t1
        print('Accumulated CPU_TIME = %.2f seconds = %s (hh:mm:ss)' \
              % (cpu_time, datetime.timedelta(seconds=cpu_time)))

        count += 1

    return pathlines

def spread_pathlines(fill_args, flow, factor, printIt):
    # fill_args = births, lifetime, initial_positions
    pathlines = fill_pathlines(*fill_args)

    return compute_pathlines(pathlines, flow, factor, printIt)