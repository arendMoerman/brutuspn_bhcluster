#include <iostream>
using namespace std;
#include <vector>
#include <cmath>

#ifndef __Products_h
#define __Products_h

class Products {
  double Resc;
  double fiso_bin, fiso_flyby;

  int Nesc;
  vector<int> index_esc;
  vector< vector<double> > data_esc; 

  // 2-body binary data
  int Nbin;
  vector< vector<int> > index_bin; 
  vector< vector<double> > data_bin;
  int Nbin_esc;
  vector< vector<int> > index_bin_esc;
  vector< vector<double> > data_bin_esc;

  // 2-body flyby data
  int Nflyby;
  vector< vector<int> > index_flyby; 
  vector< vector<double> > data_flyby;  

  int Ntri_esc;
  vector< vector<int> > index_tri_esc;  
  vector< vector<double> > data_tri_esc;

  bool isDissolved;
  int outcome;

  public:

  Products();
  Products(double fiso);

  vector<double> get_cm_multiple(vector<double> &data, vector<int> &index);
  vector<double> get_coordinates_multiple(vector<double> &data, vector<int> &index);
  vector<double> get_coordinates_multiple_in_cm_frame(vector<double> &data, vector<int> &index);

  double get_kinetic_energy_particle(vector<double> &data, int index);
  double get_potential_energy_particle(vector<double> &data, int index);
  double get_energy_particle(vector<double> &data, int index);
  double get_external_kinetic_energy_multiple(vector<double> &data, vector<int> index);
  double get_external_potential_energy_multiple(vector<double> &data, vector<int> index);
  double get_external_energy_multiple(vector<double> &data, vector<int> index);
  double get_internal_kinetic_energy_multiple(vector<double> &data, vector<int> index);
  double get_internal_potential_energy_multiple(vector<double> &data, vector<int> index);
  double get_internal_energy_multiple(vector<double> &data, vector<int> index);

  double get_angular_momentum_particle(vector<double> &data, int index);
  double get_external_angular_momentum_multiple(vector<double> &data, vector<int> &index);
  double get_internal_angular_momentum_multiple(vector<double> &data, vector<int> &index);

  bool isIsolated_abs_particle(vector<double> &data, int index, double R);
  bool isIsolated_abs_multiple(vector<double> &data, vector<int> &index, double R);
  bool isIsolated_rel_multiple(vector<double> &data, vector<int> &index, double f);

  bool is_nn_multiple(vector<double> &data, vector<int> index);

  vector<int> get_nn_list_particle(vector<double> &data, int index, int numNN);
  vector<int> get_nn_list_multiple(vector<double> &data, vector<int> index, int numNN);
  double get_distance_nn_particle(vector<double> &data, int index);
  double get_distance_nn_multiple(vector<double> &data, vector<int> index);

  vector<double> get_escaper_properties(vector<double> &data, int index1);
  vector<double> get_binary_properties(vector<double> &data, int index1, int index2);
  vector<double> get_triple_properties(vector<double> &data, int index1, int index2, int index3);

  vector<double> get_binary_internal_properties(vector<double> &data, int index1, int index2);

  //-----------------------------------------------

  // Detect isolated, escapers
  void detect_escapers(vector<double> &data, double Resc);

  // Detect isolated binaries
  void detect_binaries(vector<double> &data);
  void detect_binaries(vector<double> &data, double Resc, double fiso_bin);

  // Detect isolated 2-body encounters
  void detect_2body_encounters(vector<double> &data);

  // Detect isolated 3-body encounters
  void detect_3body_encounters(vector<double> &data);

  // Detect isolated 4-body encounters
  void detect_4body_encounters(vector<double> &data);

  //-----------------------------------------------

  void process_4body(vector<double> &data);

  //-----------------------------------------------

  int get_Nesc();
  vector<int> get_index_esc();
  vector< vector<double> > get_data_esc(); 

  int get_Nbin();
  vector< vector<int> > get_index_bin();
  vector< vector<double> > get_data_bin(); 

  int get_Nbin_esc();
  vector< vector<int> > get_index_bin_esc();
  vector< vector<double> > get_data_bin_esc(); 

  int get_Nflyby();
  vector< vector<int> > get_index_flyby();
  vector< vector<double> > get_data_flyby(); 

  int get_Ntri_esc();
  vector< vector<int> > get_index_tri_esc();
  vector< vector<double> > get_data_tri_esc(); 

  bool get_isDissolved();
  int get_outcome();

  void remove_particles(vector<double> &data, vector<int> &index);
  void add_particles(vector<double> &data, vector<double> &d);
  void add_particles(vector<double> &data, vector< vector<double> > &d);
  void replace_cm_particles(vector<double> &data, vector< vector<double> > &coordinates_bin);
};

#endif


