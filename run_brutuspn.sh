#!/bin/bash

# USER PARAMETERS
# Name of initial condition
INITCOND="Nbh_10_run_1"

# Start time in kyr
STARTTIME=0

# Endtime in kyr
ENDTTIME=5

# Snapshot interval in kyr
DTSNAP=0.00001

# Positive power of tolerance of Bulirsch-Stoer integration, converted internally to 10**-(power)
EPS=10

################################################################################################
# The following parameters usually do not have to be fine-tuned, they are ususally good as-is
# BrutusPN executable
EXEC="./brutuspn/main.exe"

# Number of bytes for mpreal, given EPS
LW=$((4 * EPS + 32)) 

# Initial Bulirsch-Stoer timestep in nbody units
ETA=0.1

# Max number of Bulirsch-Stoer steps before non-convergence
NMAX=64

# Location of initial condition to run
INITLOC="./brutuspn/inputs/"${INITCOND}".txt.in"

# Output location of simulation. Automatically placed in ./outputs/, default is same name is init cond.
OUTFILE=${INITCOND}

# Read first line of inital condition to obtain number of bodies in simulation
LINE0=($(head -n 1 ${INITLOC}))

# Number of bodies in simulation
NBODIES=${LINE0[1]}

# Mode for inital condition (for the cluster simulations, leave this at 'file')
MODE="file"
echo "Running initial condition:" ${INITLOC} ${NBODIES}

${EXEC} ${OUTFILE} ${STARTTIME} ${ENDTTIME} ${DTSNAP} ${ETA} ${EPS} ${LW} ${NMAX} ${NBODIES} ${MODE} ${INITLOC} 

