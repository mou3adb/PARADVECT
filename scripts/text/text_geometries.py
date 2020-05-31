"""
The geometries file has the following structure:

N_geometries
(blank line)
N_nodes_geometry1
node0
node1
...
(blank line)
N_nodes_geometry2
node0
node1
...
(blank line)

...
...
...

N_nodes_geometryG
node0
node1
...

"""
import numpy as np
#==============================================================================
def write_geometries(geometries_nodes, outfile):
    f = open(outfile, 'w')

    print('Writing geometries...')

    N_geometries = len(geometries_nodes)
    f.write('%d\n' % N_geometries)

    f.write('\n') # blank line

    for i in range(N_geometries):
        geometry_nodes = geometries_nodes[i]

        N_nodes_of_geometry_i = len(geometry_nodes)

        f.write('%d\n' % N_nodes_of_geometry_i)

        str_nodes = '\n'.join([str(g) for g in geometry_nodes]) + '\n'
        f.write(str_nodes)

        f.write('\n')

    f.close()

def read_geometries(infile):
    f = open(infile, 'r')

    geometries_nodes = []

    N_geometries = int(f.readline().strip())

    f.readline() # blank line

    for i in range(N_geometries):
        N_nodes = int(f.readline().strip())

        geometry_nodes = []

        for j in range(N_nodes):
            geometry_nodes.append(int(f.readline().strip()))

        geometries_nodes.append(geometry_nodes)

        f.readline() # blank line

    f.close()

    return np.array(geometries_nodes)

