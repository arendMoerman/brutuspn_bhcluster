import numpy as np
import argparse

MSOL = 1.9884e30  #kg
PARSEC = 3.086e16 #meters
G = 6.67430e-11 #m**3 / kg / (s**2)

VEL = np.sqrt(G * MSOL / PARSEC) # m/s
ZETA = 299792458 / VEL
TIME = PARSEC / VEL / 60 / 60 / 24 / 365 / 1000 # kyr

def get_zeta_timescale():
    parser = argparse.ArgumentParser(description="options for calculating zeta and timescaling from mass and distance scale")
    parser.add_argument("-m", "--mass", type=float, help="Set mass scale in solar masses")
    parser.add_argument("-r", "--distance", type=float, help="Set distance scale in parsecs")
    args = parser.parse_args()

    fm = parser.mass
    print(fm)

if __name__ == "__main__":
    get_zeta_timescale()
