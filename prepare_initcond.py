import numpy as np
import csv

MsgrA = 4.297e6 #Msol, Mass scale (Mass from https://www.aanda.org/articles/aa/full_html/2023/09/aa47416-23/aa47416-23.html)
RNUM = 0.5 # Distance scale in parsecs

MSOL = 1.9884e30  #kg
PARSEC = 3.086e16 #meters
G = 6.67430e-11 #m**3 / kg / (s**2)

FM = MSOL * MsgrA
FR = PARSEC * RNUM 

VEL = np.sqrt(G * FM / FR) # m/s
ZETA = 299792458 / VEL
TIME = FR / VEL / 60 / 60 / 24 / 365 / 1000 # kyr

def convert(name):
    Nbodies = 0
    with open(f'Initial_Conditions/{name}', newline='')  as csvfile, open(f"brutuspn/inputs/{name}.in", "w") as outfile:
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
