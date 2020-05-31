import time, datetime
import numpy as np

from pathline import Pathline
#==============================================================================
class Frontline(object):
    """
    Frontline class calculate the frontline of a flow

    Attributes
    ----------

    segment: numpy array of floats
        The ordinates from which you launch your hydrogen bubbles or dye.

    birth: integer
        The time when you start to launch your dye. It is relative to the
        timeline of the flow.
        e.g. birth = 0 is the beginning of the flow.
             birth = int(len(flow.times)/2.) is the middle of the flow

    lifetime: integer
        The lifetime of your dye.
        Should be in the same units of the flow timeline.

    front: numpy array of pathline
        The frontline is by definition the first set of bubble that are launched
        at the same time. These are basically pathlines started at the same time.

    """
    def __init__(self, segment, birth, lifetime):
        self.segment = segment

        self.birth = birth
        self.lifetime = lifetime

        self.front = []

    def compute_frontline(self, flow, factor, printIt):
        t1 = time.time()

        for pos0 in self.segment:
            pathline = Pathline(self.birth, self.lifetime, pos0)

            pathline.compute_pathline(flow, factor, printIt)

            self.front.append(pathline.path)

        self.front = np.array(self.front)

        d = time.time() - t1
        print('===/ Frontline successfully computed :)')
        print('===/ CPU_TIME = %.2f seconds = %s (hh:mm:ss)\n' \
              % (d, datetime.timedelta(seconds=d)))