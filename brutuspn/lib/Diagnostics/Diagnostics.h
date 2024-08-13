using namespace std;
#include <vector>
#include <cmath>

#include "Tools.h"

#ifndef __Diagnostics_h
#define __Diagnostics_h

class Diagnostics {
  public:

  // Conserved quantities
  double get_mass(vector<double> &data);

  vector<double> get_rcm(vector<double> &data);
  vector<double> get_vcm(vector<double> &data);
  vector<double> get_lcm(vector<double> &data);

  double get_kinetic_energy(vector<double> &data);
  double get_potential_energy(vector<double> &data);
  double get_energy(vector<double> &data);

  // System properties
  double get_virial_radius(vector<double> &data);
  double get_harmonic_radius(vector<double> &data);
  double get_velocity_disperion(vector<double> &data);

  // Core collapse properties
  vector<double> get_rho(vector<double> &data);
  vector<double> get_rho_center(vector<double> &data, vector<double> &rho);
  double get_core_radius85(vector<double> &data, vector<double> &rho, vector<double> &rho_center);
  double get_core_radius90(vector<double> &data, vector<double> &rho, vector<double> &rho_center);  
  double get_core_density(vector<double> &data, vector<double> &rho);
  double get_core_mass(vector<double> &data, double &core_radius);
  int get_core_number(vector<double> &data, double &core_radius);
};

#endif

