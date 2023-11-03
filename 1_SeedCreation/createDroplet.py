 
#!/usr/bin/env python3

import numpy as np
import glob
import os
import os.path
import sys

import ovito
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *

# python3 createDroplet.py droplet.data Droplet 500
cavityFile = sys.argv[1]
outDir = sys.argv[2]
numAtoms = int(sys.argv[3])

if not numAtoms % 2 == 0:
    print('Number of atoms to remove is not even')
    exit()

########################  Create droplet  ################################
# Do read database
output_file = 'droplet'+str(numAtoms)+'.data'
output_file = outDir+'/'+str(numAtoms)+'/'+output_file
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Compute distance to zero position
pipeline = import_file(cavityFile, multiple_frames=True)
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='distanceCenter',
    expressions=['sqrt((Position.X)^2 + (Position.Y)^2 + (Position.Z)^2)']
))
data = pipeline.compute(0)
#Create list with id's, distances and atom type
distances = []
for i in range(len(data.particles.get('distanceCenter'))):
    distances.append([data.particles.get('distanceCenter')[i], data.particles.get('Particle Identifier')[i], data.particles.get('Particle Type')[i]])
#Sort by distance
sortedDistances = sorted(distances)
print(sortedDistances[0])
#Select atoms closest to center
deleteZn = []
deleteO = []
for i in range(len(sortedDistances)):
    if sortedDistances[i][2] == 1:
        if len(deleteO) < numAtoms / 2:
            deleteO.append(sortedDistances[i][1])
    if sortedDistances[i][2] == 2:
        if len(deleteZn) < numAtoms / 2:
            deleteZn.append(sortedDistances[i][1])
atomsDelete = deleteZn + deleteO
#Delete surrounding
ids = data.particles["Particle Identifier"]
data.particles_.delete_elements(np.in1d(ids, atomsDelete, assume_unique = True, invert = True))
#Export data file
export_file(data, output_file, "lammps/data", atom_style="charge", ignore_identifiers=True, frame=pipeline.source.num_frames)