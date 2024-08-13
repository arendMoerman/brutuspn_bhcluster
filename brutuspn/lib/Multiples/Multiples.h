#include <vector>
#include <cmath>

#ifndef __MULTIPLES_H
#define __MULTIPLES_H

class Multiples {
  double R_encounter;
  double f_isolated;

  vector<int> nn_index;
  vector<double> nn_distance;
  vector<double> nn2_distance;

  vector< vector<int> > two_body_flybys;
  vector< vector<int> > two_body_binaries;

  public:

  Multiples();

  // 2-body encounters
  void make_nn_list(vector<double> &data);
  vector<double> get_binary_properties(vector<double> &data, int index1, int index2);
  void detect_two_body_encounters(vector<double> &data);

  // 3-body encounters


  // 4-body encounters

  // total
  void process(vector<double> &data);

};

#endif


