#!/bin/bash

echo "Submitting jobs for different random seeds"

export randomSeed=00000
sbatch xminim

export randomSeed=11111
sbatch xminim

export randomSeed=22222
sbatch xminim

export randomSeed=33333
sbatch xminim

export randomSeed=44444
sbatch xminim
