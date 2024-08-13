/*
Potential types:
- none (default)
- point_mass
*/

#include <iostream>
using namespace std;

#include <cmath>
#include <string>
#include <vector>

#ifndef __Potentials_h
#define __Potentials_h

class Potentials {
  string potential_type;
  double tolerance;

  double Mgal;
  vector<double> par;

  public:

  Potentials();

  void set_potential_type(string type, vector<double> parm);
  void set_tolerance(double tol);

  double get_point_mass_potential_at(double m, double x, double y, double z);
  double get_homo_sphere_potential_at(double m, double x, double y, double z);
  double get_isochrone_potential_at(double m, double x, double y, double z);

  vector<double> get_acceleration_at(double m, double x, double y, double z);

  vector<double> iterate(double phi0, double dq, double m, double x, double y, double z);
  void extrapol(vector<double> &dqs, vector< vector<double> > &as, vector<double> &a_exp);
  double extrapolate(vector<double> x, vector<double> y, double x0);
  bool error_control(vector<double> &a1, vector<double> &a2);
};

#endif


