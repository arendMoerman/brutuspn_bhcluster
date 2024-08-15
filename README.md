# Installation of BrutusPN
## Prerequisites
- Arbitrary-precision libraries:

  gmp		: https://gmplib.org/
  mpfr		: http://www.mpfr.org/
  mpfr c++	: http://www.holoborodko.com/pavel/mpfr/

- Build the code

  1) Adjust the path to brutuspn_bhcluster/brutuspn/mpreal/mpreal.h in makefile
  2) Type make

# Instructions for starting a BrutusPN simulation for the BH cluster experiment

## Step 1. set mass and distance scale
The mass and distance scale are set by adjusting the ```MNUM``` and ```RNUM``` parameters in the ```prepare_initcond.py``` file, respectively.
The numbers are fractions of the mass of SgrA* [Msol] and a parsec, respectively.
After setting, run the python script and all initial conditions will be converted to the Nbody units.
Also, the dimensionless speed of light will be set and the internal timescaling will be set.

