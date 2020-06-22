"""
In this file, compute_args = flow, factor, printIt

SEMANTIC NOTE: 'precision' should be replaced by 'resolution'.
  In fact, precision means repeatability, whereas resolution is the smallest
  increment. Accuracy is another word, but means how close a measure is to
  the true value. Error is not the right term either, because it tells the
  difference between the measured and true value. In our case, the
  uncertainty is linked to the resolution, e.g. plus or minus diameter/1000.

"""
import time, datetime

import numpy as np

from copy import deepcopy

from particle import Particle
#==============================================================================
def dichotomy(p1, p2, precision, flow, factor, printIt):
    """
    This function returns the ordinate of the capture region boundary.

    p1 and p2 should already be computed, and one of the particles should
    be a captured one.
    """
    x0 = p1.pos0[0]
    y1, y2 = p1.pos0[1], p2.pos0[1]

    epsilon = precision*p1.diameter

    # Check whether p1 and p2 have different birth dates and lifetimes
    if p1.birth != p2.birth:
        exception_msg = 'The two particles have different birth dates (-_-). '\
                      + 'They need to have the same.'
        raise Exception(exception_msg)

    if p1.lifetime != p2.lifetime:
        exception_msg = 'The two particles have different lifetimes (-_-). '\
                      + 'They need to have the same.'
        raise Exception(exception_msg)

    if p1.captured == p2.captured:
        if p1.captured:
            exception_msg = 'Both are captured (-_-). Please try other initial '\
                          + 'positions.'
            raise Exception(exception_msg)

        else:
            exception_msg = 'Both are uncaptured (-_-). Please try other initial '\
                          + 'positions.'
            raise Exception(exception_msg)

    elif abs(y2 - y1) < epsilon:
        print('Resolution: epsilon = precision*R = %e reached' % epsilon)

        if p1.captured:
            print('y_limit = %e' % y1)
            return y1

        else:
            print('y_limit = %e' % y2)
            return y2

    else:
        print('Resolution not reached yet')
        print('Dichotomy between y = %e and y = %e' % (min(y1,y2),max(y1,y2)))

        ymean = 0.5*(y1 + y2)

        new_pos0 = np.array([x0, ymean])

        pmean = Particle(p1.diameter, p1.density, p1.birth, p1.lifetime,
                         new_pos0, None)

        pmean.compute_trajectory(flow, factor, printIt)

        if pmean.captured:
            if not(p1.captured):
                return dichotomy(p1, pmean, precision, flow, factor, printIt)
            if not(p2.captured):
                return dichotomy(p2, pmean, precision, flow, factor, printIt)

        else:
            if p1.captured:
                return dichotomy(p1, pmean, precision, flow, factor, printIt)
            if p2.captured:
                return dichotomy(p2, pmean, precision, flow, factor, printIt)

def particle_up(p):
    k = 0
    flag, bary_pos = p.is_up[0]

    # Increment k until the particle attains the object, ie. until flag != None.
    # The second condition is when the particle can't even approach the object.
    while flag == None and k < len(p.is_up) - 1:
        k += 1
        flag, bary_pos = p.is_up[k]

    if flag == None:
        exception_msg = 'The particle hasn\'t reached the object. ' \
                      + 'Please try to increase the lifetime of the particle.'
        raise Exception(exception_msg)

    flag1, bary_pos1 = flag, bary_pos

    # In this case, we are sure flag != None. The particle is either above or
    # below the object.
    while flag != None and k < len(p.is_up) - 1:
        k += 1
        flag, bary_pos = p.is_up[k]

    if k == len(p.is_up) - 1:
        # This means that the particle is still in the same upper or
        # lower part of the object.
        return flag1

    else:
        assert(flag == None)
        if bary_pos*bary_pos1 < 0:
            # It this case the particle has moved from above to
            # below the object.
            return flag1

        else:
            # In this case the particle has returned back.
            #
            # WARNING:
            # This assumes that the particle would not return back a
            # second time and change its emplacement
            # e.g.
            # If the particle was going upper the object, then returned
            # back in the front and started to go underneath, we assume that it
            # continues to go underneath and doesn't return back again to
            # switch to the upper part.
            pCopy = deepcopy(p)
            pCopy.is_up = pCopy.is_up[k:]

            return particle_up(pCopy)

def capture_ordinate(p, increment, flow, factor, printIt):
    """
    This function finds a captured particles.

    It take a particle p, then sees:
        if p is captured, that's the particle we were looking for,
        if not, launch a particle just underneath and see again.
    """
    x0, y0 = p.pos0

    forbid_up, forbid_down = 1., -1.
    if not(forbid_down < y0 < forbid_up):
        exception_msg = 'The starting point of the particle is outside the ' \
                      + 'prescribed limits. Please try to start close to ' \
                      + 'the capture domain.'
        raise Exception(exception_msg)

    if p.captured:
        return p, increment

    else:
        print('Calculating the capture ordinate, starting from y0 = %e' % y0)

        is_up = particle_up(p)

        if (increment > 0 and is_up) or (increment < 0 and not(is_up)):

            print('Still going above the object... Repeat again.')

            new_pos0 = np.array([x0, y0 - increment])

            pNext = Particle(p.diameter, p.density, p.birth, p.lifetime,
                             new_pos0, None)

            pNext.compute_trajectory(flow, factor, printIt)

            return capture_ordinate(pNext, increment, flow, factor, printIt)

        else:
            print('Now it\'s going below the object!')

            new_pos0 = np.array([x0, y0 + 0.7*increment])

            pNext = Particle(p.diameter, p.density, p.birth, p.lifetime,
                             new_pos0, None)

            pNext.compute_trajectory(flow, factor, printIt)

            return capture_ordinate(pNext, 0.3*increment, flow, factor, printIt)

def search_the_limit_outwards(p, increment, flow, factor, printIt):
    """
    p should be captured.

    If increment > 0, it searches the lower limit.
    """
    if not(p.captured):
        exception_msg = 'Sorry, the particle isn\'t captured. We cannot ' \
                      + 'proceed to the computation :('
        raise Exception(exception_msg)

    x0, y0 = p.pos0

    new_p = deepcopy(p)

    while new_p.captured:
        print('y0 = %e. Still inside the capture domain...' % new_p.pos0[1])

        pCopy = deepcopy(new_p)

        y0 -= increment
        new_pos0 = np.array([x0, y0])

        new_p = Particle(p.diameter, p.density, p.birth, p.lifetime,
                         new_pos0, None)

        new_p.compute_trajectory(flow, factor, printIt)

    return pCopy, new_p

def find_limits(p, increment, precision, flow, factor, printIt):
    if p.captured:
        pCaptured, finalIncrement = p, increment

    elif particle_up(p):
        print('Particle goes above the object.')
        pCaptured, finalIncrement = capture_ordinate(p, increment,
                                                     flow, factor, printIt)

    else:
        print('Particle goes below the object. Everything is reversed upside down.')
        pCaptured, finalIncrement = capture_ordinate(p, -increment,
                                                     flow, factor, printIt)

    yCapture = pCaptured.pos0[1]
    print('')
    print(40*'*')
    print('GOTCHA!')
    print('y_capture       = %e' % yCapture)
    print('Final increment = %e' % finalIncrement)
    print(40*'*')

    print('')
    print(40*'*')
    print('Searching the limit upwards')
    print(40*'*')

    pCapturedUp, pUp = search_the_limit_outwards(pCaptured, -0.3*finalIncrement,
                                                 flow, factor, printIt)

    print('')
    print(40*'*')
    print('REMINDER')
    print('%s, R = %s' % (flow, p.diameter))
    print(40*'*')

    print('')
    print(40*'*')
    print('Searching the limit downwards')
    print(40*'*')
    pCapturedDown, pDown = search_the_limit_outwards(pCaptured, 0.3*finalIncrement,
                                                     flow, factor, printIt)

    print('')
    print(40*'*')
    print('REMINDER')
    print('%s, R = %s' % (flow, p.diameter))
    print(40*'*')

    print('')
    print(40*'*')
    print('Starting dichotomy for the upper bound')
    print(40*'*')
    limitUp = dichotomy(pCaptured, pUp, precision, flow, factor, printIt)

    print('')
    print(40*'*')
    print('Starting dichotomy for the lower bound')
    print(40*'*')
    limitDown = dichotomy(pDown, pCaptured, precision, flow, factor, printIt)

    return limitUp, limitDown

def automatic_search(p, increment, precision, flow, factor, printIt):
    t1 = time.clock()

    limUp, limDown = find_limits(p, increment, precision, flow, factor, printIt)

    e = limUp - limDown

    d = time.clock() - t1
    print('')
    print('FINALLY, WE GOT THE LIMITS! :)')
    print('CPU_TIME = %.2f seconds = %s (hh:mm:ss)' % (d, datetime.timedelta(seconds=d)))
    print('')

    return e, (limDown, limUp)
