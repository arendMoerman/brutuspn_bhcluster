# Installation of BrutusPN
## Prerequisites
- Arbitrary-precision libraries:

  gmp		: https://gmplib.org/
  mpfr		: http://www.mpfr.org/
  mpfr c++	: http://www.holoborodko.com/pavel/mpfr/

- Build the code

  1) Adjust the path to ```brutuspn_bhcluster/brutuspn/mpreal/mpreal.h``` in makefile
  2) Change working directory to ```/brutuspn/``` directory
  3) Type ```make```

# Instructions for starting a BrutusPN simulation for the BH cluster experiment

## Step 0. set mass and distance scale, run ```prepare_initcond.py```
The mass and distance scale are set by adjusting the ```MNUM``` and ```RNUM``` parameters in the ```prepare_initcond.py``` file, respectively.
The numbers are fractions of the mass of SgrA* [Msol] and a parsec, respectively.
After setting, run the python script and all initial conditions will be converted to the Nbody units.
Also, the dimensionless speed of light will be set and the internal timescaling will be set.
Remember, this script must be run at least ONCE on every new machine that the simulations will be run on.

## Step 1. Adjust parameters in ```run_brutuspn.sh```
Open the ```run_brutuspn.sh``` file and adjust the following to your requirements:
* ```INITCOND```: Name of initial condition to run (without ANY extensions)
* ```STARTTIME```: Starting time of simulation in kyr (can usually be kept at zero)
* ```ENDTTIME```: Endtime of simulation in kyr
* ```DTSNAP```: Snapshot interval for writing output in kyr
* ```EPS```: Power of Bulirsch-Stoer tolerance. Note that it should be positive. For example, if a tolerance of 1e-10 is required, set ```EPS=10```.
The other parameters can be changed but usually dont need changing.

## Step 2. Run ```bash run_brutuspn.sh```
The output will be stored in a folder called ```./outputs/```, in the following format:
* ```<INITCOND>.out```: Position and velocities of each body, as function of time
* ```<INITCOND>.energies```: Total energy of the system as function of time. This file contains the timestamps
* ```<INITCOND>.log```: This file contains a log of the simulation. TODO: include keys of mergers in this file.

## Step 3. Inspect the output using ```plot_output.py```
The output can be inspected in real-time (i.e., during a simulation), by running the following command: ```python plot_output.py -n <INITCOND>```.
The script automatically plots the simulation projected in the x-y plane, but this can be changed by passing the extra ```-p``` command. 
For an overview of projection options, run ```python plot_output.py -h```.

## TODO
Add a parameter to set merging radius. Currently, ```BrutusPN``` assumes a merger within one R_ISCO, but I want to add a parameter, call it ```alpha```, to set the merger radius to ```alpha```*R_ISCO.
Also, need to add merger keys to the .log file.
