#include <iostream>
using namespace std;
#include <vector>
#include <cmath>

#ifndef __Delta_h
#define __Delta_h

class Delta {
  public:

  double get_phase_space_distance(vector<double> &data1, vector<double> &data2);
  double get_phase_space_distance_normalized(vector<double> &data1, vector<double> &data2);

  double get_position_space_distance(vector<double> &data1, vector<double> &data2);
  double get_position_space_distance_normalized(vector<double> &data1, vector<double> &data2);

  double get_velocity_space_distance(vector<double> &data1, vector<double> &data2);
  double get_velocity_space_distance_normalized(vector<double> &data1, vector<double> &data2);
};

#endif


