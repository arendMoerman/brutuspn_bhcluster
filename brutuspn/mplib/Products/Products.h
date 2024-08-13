#include <iostream>
using namespace std;
#include <vector>
#include <cmath>

#include "/usr/include/mpreal.h"
using namespace mpfr;

#ifndef __Products_h
#define __Products_h

class Products {
  mpreal Resc;
  mpreal fiso_bin, fiso_flyby;

  int Nesc;
  vector<int> index_esc;
  vector< vector<mpreal> > data_esc; 

  // 2-body binary data
  int Nbin;
  vector< vector<int> > index_bin; 
  vector< vector<mpreal> > data_bin;
  int Nbin_esc;
  vector< vector<int> > index_bin_esc;
  vector< vector<mpreal> > data_bin_esc;

  // 2-body flyby data
  int Nflyby;
  vector< vector<int> > index_flyby; 
  vector< vector<mpreal> > data_flyby;  

  int Ntri_esc;
  vector< vector<int> > index_tri_esc;  
  vector< vector<mpreal> > data_tri_esc;

  bool isDissolved;
  int outcome;

  public:

  Products();
  Products(mpreal fiso);

  vector<mpreal> get_cm_multiple(vector<mpreal> &data, vector<int> &index);
  vector<mpreal> get_coordinates_multiple(vector<mpreal> &data, vector<int> &index);
  vector<mpreal> get_coordinates_multiple_in_cm_frame(vector<mpreal> &data, vector<int> &index);

  mpreal get_kinetic_energy_particle(vector<mpreal> &data, int index);
  mpreal get_potential_energy_particle(vector<mpreal> &data, int index);
  mpreal get_energy_particle(vector<mpreal> &data, int index);
  mpreal get_external_kinetic_energy_multiple(vector<mpreal> &data, vector<int> index);
  mpreal get_external_potential_energy_multiple(vector<mpreal> &data, vector<int> index);
  mpreal get_external_energy_multiple(vector<mpreal> &data, vector<int> index);
  mpreal get_internal_kinetic_energy_multiple(vector<mpreal> &data, vector<int> index);
  mpreal get_internal_potential_energy_multiple(vector<mpreal> &data, vector<int> index);
  mpreal get_internal_energy_multiple(vector<mpreal> &data, vector<int> index);

  mpreal get_angular_momentum_particle(vector<mpreal> &data, int index);
  mpreal get_external_angular_momentum_multiple(vector<mpreal> &data, vector<int> &index);
  mpreal get_internal_angular_momentum_multiple(vector<mpreal> &data, vector<int> &index);

  bool isIsolated_abs_particle(vector<mpreal> &data, int index, mpreal R);
  bool isIsolated_abs_multiple(vector<mpreal> &data, vector<int> &index, mpreal R);
  bool isIsolated_rel_multiple(vector<mpreal> &data, vector<int> &index, mpreal f);

  bool is_nn_multiple(vector<mpreal> &data, vector<int> index);

  vector<int> get_nn_list_particle(vector<mpreal> &data, int index, int numNN);
  vector<int> get_nn_list_multiple(vector<mpreal> &data, vector<int> index, int numNN);
  mpreal get_distance_nn_particle(vector<mpreal> &data, int index);
  mpreal get_distance_nn_multiple(vector<mpreal> &data, vector<int> index);

  vector<mpreal> get_escaper_properties(vector<mpreal> &data, int index1);
  vector<mpreal> get_binary_properties(vector<mpreal> &data, int index1, int index2);
  vector<mpreal> get_triple_properties(vector<mpreal> &data, int index1, int index2, int index3);

  vector<mpreal> get_binary_internal_properties(vector<mpreal> &data, int index1, int index2);

  //-----------------------------------------------

  // Detect isolated, escapers
  void detect_escapers(vector<mpreal> &data, mpreal Resc);

  // Detect isolated binaries
  void detect_binaries(vector<mpreal> &data);
  void detect_binaries(vector<mpreal> &data, mpreal Resc, mpreal fiso_bin);

  // Detect isolated 2-body encounters
  void detect_2body_encounters(vector<mpreal> &data);

  // Detect isolated 3-body encounters
  void detect_3body_encounters(vector<mpreal> &data);

  // Detect isolated 4-body encounters
  void detect_4body_encounters(vector<mpreal> &data);

  //-----------------------------------------------

  void process_4body(vector<mpreal> &data);

  //-----------------------------------------------

  int get_Nesc();
  vector<int> get_index_esc();
  vector< vector<mpreal> > get_data_esc(); 

  int get_Nbin();
  vector< vector<int> > get_index_bin();
  vector< vector<mpreal> > get_data_bin(); 

  int get_Nbin_esc();
  vector< vector<int> > get_index_bin_esc();
  vector< vector<mpreal> > get_data_bin_esc(); 

  int get_Nflyby();
  vector< vector<int> > get_index_flyby();
  vector< vector<mpreal> > get_data_flyby(); 

  int get_Ntri_esc();
  vector< vector<int> > get_index_tri_esc();
  vector< vector<mpreal> > get_data_tri_esc(); 

  bool get_isDissolved();
  int get_outcome();

  void remove_particles(vector<mpreal> &data, vector<int> &index);
  void add_particles(vector<mpreal> &data, vector<mpreal> &d);
  void add_particles(vector<mpreal> &data, vector< vector<mpreal> > &d);
  void replace_cm_particles(vector<mpreal> &data, vector< vector<mpreal> > &coordinates_bin);
};

#endif


