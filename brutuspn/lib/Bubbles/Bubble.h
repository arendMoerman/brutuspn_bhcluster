using namespace std;
#include <vector>
#include <cmath>

#ifndef __BUBBLE_H
#define __BUBBLE_H

class Bubble {
  int N;			// Number of stars
  vector<double> data;		// Coordinates relative to CoM

  vector<double> data_cm;	// Coordinates CoM particle

  double R;			// Bubble radius
  double f_radius;		// Bubble radius parameter

  double dt;			// Timescale

  public:

  Bubble();
  Bubble(double &f_radius);
  Bubble(vector<double> &raw_data);  
  Bubble(double &f_radius, vector<double> &raw_data);  

  void set_N(int &N);
  void set_data(vector<double> &data); 
  void set_data_cm(vector<double> &data); 
  void set_R(double &R);
  void set_f_radius(double &f_radius);
  void set_dt(double &dt);

  int get_N();
  vector<double> get_data();
  vector<double> get_data_cm();
  double get_R();
  double get_f_radius();
  double get_dt();

  void process_raw_data(vector<double> &raw_data);
  vector<double> get_raw_data();
};

#endif


