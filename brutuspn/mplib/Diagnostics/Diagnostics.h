using namespace std;
#include <vector>
#include <cmath>

#include "Tools.h"

#include "/usr/include/mpreal.h"
using namespace mpfr;

#ifndef __Diagnostics_h
#define __Diagnostics_h

class Diagnostics {
  public:

  // Conserved quantities
  mpreal get_mass(vector<mpreal> &data);

  vector<mpreal> get_rcm(vector<mpreal> &data);
  vector<mpreal> get_vcm(vector<mpreal> &data);
  vector<mpreal> get_lcm(vector<mpreal> &data);

  mpreal get_kinetic_energy(vector<mpreal> &data);
  mpreal get_potential_energy(vector<mpreal> &data);
  mpreal get_energy(vector<mpreal> &data);

  // System properties
  mpreal get_virial_radius(vector<mpreal> &data);
  mpreal get_harmonic_radius(vector<mpreal> &data);
  mpreal get_velocity_dispersion(vector<mpreal> &data);
  vector<mpreal> get_velocity_dispersion_1d(vector<mpreal> &data);

  // Individual quantities
  mpreal get_kinetic_energy(vector<mpreal> &data, int index);
  mpreal get_potential_energy(vector<mpreal> &data, int index);
  mpreal get_energy(vector<mpreal> &data, int index);

  mpreal get_radius(vector<mpreal> &data, int index);
  mpreal get_speed(vector<mpreal> &data, int index);

  // Core collapse properties
  vector<mpreal> get_rho(vector<mpreal> &data);
  vector<mpreal> get_rho_center(vector<mpreal> &data, vector<mpreal> &rho);
  mpreal get_core_radius85(vector<mpreal> &data, vector<mpreal> &rho, vector<mpreal> &rho_center);
  mpreal get_core_radius90(vector<mpreal> &data, vector<mpreal> &rho, vector<mpreal> &rho_center);  
  mpreal get_core_density(vector<mpreal> &data, vector<mpreal> &rho);
  mpreal get_core_mass(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal &core_radius);
  int get_core_number(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal &core_radius);
  mpreal get_halfmass_radius(vector<mpreal> &data, vector<mpreal> &rho_center);

  mpreal get_velocity_dispersion_halfmass(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rh);
  vector<mpreal> get_velocity_dispersion_halfmass_1d(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rh);
  vector<mpreal> get_velocity_dispersion_plummer_radius_1d(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rp);

  vector<mpreal> get_radii(vector<mpreal> &data, vector<mpreal> &rho_center);
};

#endif

