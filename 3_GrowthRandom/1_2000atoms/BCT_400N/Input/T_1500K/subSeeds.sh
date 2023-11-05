#!/bin/bash

echo "Submitting jobs for different random seeds"
export tempG=1500

export randomSeed=11111
sbatch xrun

export randomSeed=22222
sbatch xrun

export randomSeed=33333
sbatch xrun

export randomSeed=44444
sbatch xrun

export randomSeed=55555
sbatch xrun
