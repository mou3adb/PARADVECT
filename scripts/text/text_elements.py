"""
The elements file has the following structure:

N_elements
element0_id, node0_id node1_id ..., neighbor0_id neighbor1_id ...
element1_id, node0_id node1_id ..., neighbor0_id neighbor1_id ...
...
...
...
elementN-1_id, node0_id node1_id ..., neighbor0_id neighbor1_id ...

NOTE: there are only 2 commas per line: commas separating element ids and node
      ids, and node ids and neighbors ids.

"""
import numpy as np

from element import Element
#==============================================================================
def intIt(l):
    return np.array([int(e) for e in l])

def write_elements(elements, outfile):
    f = open(outfile, 'w')

    N_elements = len(elements)
    f.write('%d\n' % N_elements)

    for i in range(N_elements):
        e = elements[i]

        str_nodes     = ('%6d '*len(e.nodes))[:-1] % tuple(e.nodes)
        str_neighbors = ('%6d '*len(e.neighbors))[:-1]  % tuple(e.neighbors)

        f.write('%6d, %s, %s\n' % (e.id, str_nodes, str_neighbors))

    f.close()

def read_elements(infile):
    f = open(infile, 'r')

    elements = []

    N_elements = int(f.readline().strip())

    for i in range(N_elements):
        fragmented = f.readline().strip().split(', ')

        e = Element(int(fragmented[0]))

        e.nodes     = intIt(fragmented[1].strip().split())
        e.neighbors = intIt(fragmented[2].strip().split())

        elements.append(e)

    f.close()

    return np.array(elements)
