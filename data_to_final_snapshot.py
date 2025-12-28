from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.lab import *

import csv
import natsort
import glob
import pickle as pkl
import os
import numpy as np


def gw_timescale(m1, m2, eccentricity, semimajor):
    """
    Function to calculate the gravitational wave timescale using Peters 1963 formula
    
    Args:
        m1 (float): Mass of the primary object
        m2 (float): Mass of the secondary object
        eccentricity (float): Eccentricity of the binary
        semimajor (float): Semi-major axis of the binary
    Returns:
        tgw (float): Gravitational wave inspiral time
    """
    mu = (m1*m2)/(m1+m2)
    coeff = (5.*constants.c**5.)/(256.*constants.G**3.)
    tgw = coeff*(semimajor**4.*(1.-eccentricity**2.)**(7./2.))/(mu*(m1+m2)**2.)
    return tgw

def store_data(new_particle, old_particle):
    """
    Transfer data onto an AMUSE particle set
    
    Args:
        new_particle (object): New particle to store data
        old_particle (object): Old particle to transfer data from
    """
    new_particle.mass = old_particle[1]
    new_particle.x = old_particle[2][0]
    new_particle.y = old_particle[2][1]
    new_particle.z = old_particle[2][2]
    new_particle.vx = old_particle[3][0]
    new_particle.vy = old_particle[3][1]
    new_particle.vz = old_particle[3][2]

folder = "rc_0.25_4e6"  # Look at MSMBH=4e6 MSun
snapshots = natsort.natsorted(glob.glob(folder+"/GRX/particle_trajectory/*"))
run_no = 0

for s in snapshots:
    with open(s, 'rb') as input_file:
        file_size = os.path.getsize(s)
        if file_size < 2.9e9:  # Larger than this, PC crashes
            print("Reading: {:}".format(input_file))
            particle_snapshot = pkl.load(input_file)
            
            final_time = particle_snapshot.iloc[:,-1]
            pre_final_time = particle_snapshot.iloc[:,-2]
            no_particles = np.shape(final_time)[0]
            particles = Particles(no_particles)
            merger_bin = Particles(2)
            
            if np.isnan(final_time.iloc[0][0]):  # Check for merger
                run_no += 1
                IC_file_name = "Final_Inspiral/Initial_Conditions/Nbh_{:}_run_{:}.txt".format(no_particles, run_no)
                merge_param_file_name = "Final_Inspiral/Merger_Parameters/Nbh_{:}_run_{:}.txt".format(no_particles, run_no)
                with open(IC_file_name, mode="w", newline="") as file:
                    brutus = [0, no_particles, 0]
                    headers = ["key", "merger [Bool]", "mass [kg]", "x [m]", "y [m]", "z [m]", "vx [m/s]", "vy [m/s]", "vz [m/s]"]
                    writer = csv.writer(file)
                    writer.writerow(brutus)
                    writer.writerow(headers)
                
                for p in range(len(final_time)):
                    particle = pre_final_time.iloc[p]

                    if np.isnan(final_time.iloc[p][0]) or final_time.iloc[p][1] > (1e6 | units.MSun):
                        merger = 1
                        if p == 0:
                            store_data(merger_bin[0], particle)
                            key1 = particle[0]
                        else:
                            store_data(merger_bin[1], particle)
                            key2 = particle[0]
                    else:
                        merger = 0        
                    store_data(particles[p], particle)
                    particles[p].merger = merger
                    particles[p].key = particle[0]
                    
                particles.move_to_center()
                
                for particle in particles:
                    data_file = [f"{particle.key:.15e}",
                                 f"{particle.merger}",
                                 f"{particle.mass.value_in(units.kg):.15e}",
                                 f"{particle.x.value_in(units.m):.15e}",
                                 f"{particle.y.value_in(units.m):.15e}",
                                 f"{particle.z.value_in(units.m):.15e}",
                                 f"{particle.vx.value_in(units.ms):.15e}",
                                 f"{particle.vy.value_in(units.ms):.15e}",
                                 f"{particle.vz.value_in(units.ms):.15e}"
                    ]
                    
                    with open(IC_file_name, mode="a", newline="") as file:
                        writer = csv.writer(file)
                        writer.writerow(data_file)
                
                kepler_elements = orbital_elements_from_binary(merger_bin, G=constants.G)
                semimajor = kepler_elements[2]
                eccentric = kepler_elements[3]
                inclinate = kepler_elements[4]
                arg_peri = kepler_elements[5]
                asc_node = kepler_elements[6]
                true_anom = kepler_elements[7]
                tgw = gw_timescale(merger_bin[0].mass,
                                   merger_bin[1].mass,
                                   eccentric,
                                   semimajor                  
                )
                distance = (merger_bin[0].position - merger_bin[1].position).lengths()
                
                IMBH_merger = merger_bin[merger_bin.mass == merger_bin.mass.min()]
                neigh_dist = (IMBH_merger.position - particles.position).lengths()
                dist_neigh = len(neigh_dist[neigh_dist > (5*distance)])
                
                lines = ["Key 1: {:}, Key 2: {:}".format(key1, key2),
                         "Semimajor: {:}".format(semimajor.in_(units.au)),
                         "Eccentricity: {:}".format(eccentric),
                         "Inclination: {:}".format(inclinate),
                         "Arg. of Periapsis: {:}".format(arg_peri),
                         "Ascending Node: {:}".format(asc_node),
                         "True Anomaly: {:}".format(true_anom),
                         "Distance: {:}".format(distance.in_(units.au)),
                         "# Particles > 2pc away from merging IMBH: {:}".format(dist_neigh),
                         "Peters tGW: {:}".format(tgw.in_(units.yr))
                ]
                
                with open(merge_param_file_name, mode="w", newline="") as file:
                    for line in lines:
                        file.write(line)
                        file.write('\n')
                