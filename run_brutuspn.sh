#!/bin/bash

#Prompt in case not added as argument
if [ -z "$1" ]
then
    read -p "Enter cluster population: " CLUSTERPOP
    read -p "Enter run index:  " RUNIDX
else
    CLUSTERPOP=$1
    RUNIDX=$2
fi

# USER PARAMETERS
# Name of initial condition
INITCOND="Nbh_${CLUSTERPOP}_run_"${RUNIDX}

# Start time in kyr
STARTTIME=0

# Endtime in kyr
ENDTTIME=5

# Snapshot interval in kyr
DTSNAP=0.0001

# Positive power of tolerance of Bulirsch-Stoer integration, converted internally to 10**-(power)
EPS=10

# Merge radius as fraction of R_ISCO
RMERGE=1

################################################################################################
# The following parameters usually do not have to be fine-tuned, they are ususally good as-is
# BrutusPN executable
EXEC="./brutuspn/main.exe"

# Number of bits for mpreal, given EPS
LW=$((4 * EPS + 32)) 

# Initial Bulirsch-Stoer timestep in nbody units
ETA=0.1

# Max number of Bulirsch-Stoer steps before non-convergence
NMAX=64

# Location of initial condition to run
INITLOC="./brutuspn/inputs/"${INITCOND}".txt.in"

if [ ! -f "$INITLOC" ]
then
	echo "Error: Initial condition doesn't exist"
	exit;
fi

# Output location of simulation. Automatically placed in ./outputs/, default is same name is init cond.
OUTFILE=${INITCOND}
OUTDIR="./outputs/${INITCOND}.energies"
if [ -f "$OUTDIR" ]
then
	echo "Error: Simulation already exists"
	exit;
fi

# Read first line of inital condition to obtain number of bodies in simulation
LINE0=($(head -n 1 ${INITLOC}))

# Number of bodies in simulation
NBODIES=${LINE0[1]}

# Mode for inital condition (for the cluster simulations, leave this at 'file')
MODE="file"
echo "Running initial condition:" ${INITLOC} ${NBODIES}

${EXEC} ${OUTFILE} ${STARTTIME} ${ENDTTIME} ${DTSNAP} ${ETA} ${EPS} ${LW} ${NMAX} ${NBODIES} ${RMERGE} ${MODE} ${INITLOC} &

done