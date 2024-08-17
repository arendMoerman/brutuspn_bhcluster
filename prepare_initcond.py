import numpy as np
import csv
import os

MsgrA = 4.297e6 #Msol, Mass scale (Mass from https://www.aanda.org/articles/aa/full_html/2023/09/aa47416-23/aa47416-23.html)
RNUM = 0.5 # Distance scale in units of parsecs
MNUM = 1 # Mass scale in units of sgra* mass [Msol]

MSOL = 1.9884e30  #kg
PARSEC = 3.086e16 #meters
G = 6.67430e-11 #m**3 / kg / (s**2)

FM = MSOL * MsgrA * MNUM
FR = PARSEC * RNUM 

VEL = np.sqrt(G * FM / FR) # m/s
ZETA = 299792458 / VEL
TIME = FR / VEL / 60 / 60 / 24 / 365 / 1000 # kyr

PNORDER = "1 1 1 1 0 0" # Order of PN terms: 1, 1xx, 2, 2.5, 3, 3.5 (3 and 3.5 PN are not well tested!). 0 means OFF, 1 means ON

def convert(name):
    Nbodies = 0

    # Write zeta to c.par
    with open("./brutuspn/c.par", "w") as cpar:
        cpar.write(str(ZETA))

    # Write time scale to t.par, so that user does not have to give endtime and snapshot interval in nbody units, but instead in kyr
    with open("./brutuspn/t.par", "w") as tpar:
        tpar.write(str(TIME))
    
    # Write used PN terms to PN.par
    with open("./brutuspn/PN.par", "w") as PNpar:
        PNpar.write(PNORDER)
    
    outdir = "./outputs/"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)    
    
    initdir = "./brutuspn/inputs/"
    if not os.path.isdir(initdir):
        os.makedirs(initdir)    

    with open(f'Initial_Conditions/{name}', newline='')  as csvfile, open(f"{initdir}{name}.in", "w") as outfile:
        initreader = csv.reader(csvfile, delimiter=',')
        for idx, row in enumerate(initreader):
            if idx == 0:
                line_w = f"{row[0]} {row[1]} {row[2]}\n"
                Nbodies = int(row[1])
                print(f"Number of bodies = {Nbodies}")
                outfile.write(line_w)
            elif idx == 1:
                continue
            else:
                row.pop(0)
                row.pop(0)
                
                Nb_mass = float(row[0]) / FM
                Nb_x = float(row[1]) / FR
                Nb_y = float(row[2]) / FR
                Nb_z = float(row[3]) / FR
                Nb_vx = float(row[4]) / VEL
                Nb_vy = float(row[5]) / VEL
                Nb_vz = float(row[6]) / VEL
                
                line_w = f"{Nb_mass} {Nb_x} {Nb_y} {Nb_z} {Nb_vx} {Nb_vy} {Nb_vz}"
                
                if idx < (Nbodies + 2):
                    line_w += "\n"
                outfile.write(line_w)


if __name__ == "__main__":
    import os
    
    for name in os.listdir("Initial_Conditions"):
        convert(name)

    print(f"Zeta = {ZETA}")
    print(f"Timescale = {TIME} kyr")
