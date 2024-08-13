#include <iostream>
using namespace std;

#include <cmath>
#include <vector>

#include "/usr/include/mpreal.h"
using namespace mpfr;

#ifndef __KEPLER_H
#define __KEPLER_H

class Kepler {
  public:

  mpreal get_eccentricity(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_semimajor_axis(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_inclination(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_right_ascension_of_ascending_node(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_argument_of_perigee(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_true_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_eccentric_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  mpreal get_mean_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);

  mpreal get_eccentricity(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_semimajor_axis(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_inclination(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_right_ascension_of_ascending_node(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_argument_of_perigee(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_true_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_eccentric_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  mpreal get_mean_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);

  mpreal get_eccentricity(vector<mpreal> &data, int index1, int index2);
  mpreal get_semimajor_axis(vector<mpreal> &data, int index1, int index2);
  mpreal get_inclination(vector<mpreal> &data, int index1, int index2);
  mpreal get_right_ascension_of_ascending_node(vector<mpreal> &data, int index1, int index2);
  mpreal get_argument_of_perigee(vector<mpreal> &data, int index1, int index2);
  mpreal get_true_anomaly(vector<mpreal> &data, int index1, int index2);
  mpreal get_eccentric_anomaly(vector<mpreal> &data, int index1, int index2);
  mpreal get_mean_anomaly(vector<mpreal> &data, int index1, int index2);

  vector<mpreal> get_orbelem(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz);
  vector<mpreal> get_orbelem(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2);
  vector<mpreal> get_orbelem(vector<mpreal> &data, int index1, int index2);

  mpreal convert_eccentric_to_true_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);
  mpreal convert_true_to_eccentric_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);
  mpreal convert_eccentric_to_mean_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);
  mpreal convert_mean_to_eccentric_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);
  mpreal convert_mean_to_true_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);
  mpreal convert_true_to_mean_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A);

  vector<mpreal> get_cartesian_coordinates(mpreal &mu, mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &E);
};

#endif
