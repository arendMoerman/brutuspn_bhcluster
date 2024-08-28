#include "Timer.h"

#include "Star.h"
#include "Cluster.h"
#include "Bulirsch_Stoer.h"

#ifndef __Brutus_h
#define __Brutus_h

class Brutus {
  mpreal t;
  int N;  
  vector<mpreal> data;

  mpreal tolerance;
  int numBits;

  mpreal eta, dt;

  Cluster cl;
  Bulirsch_Stoer bs;

  public:

  vector<mpreal> get_dt_step();
  vector<mpreal> get_tcpu_step();

  Brutus(vector<mpreal> &data, mpreal r_merge);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, mpreal r_merge);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal r_merge);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, mpreal r_merge);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, int &nmax, mpreal r_merge);

  mpreal get_eta();

  void setup(mpreal r_merge);
  void setup(int nmax, mpreal r_merge);

  void evolve(mpreal t_end, mergerOut &merge);
  void reverse_velocities();

  mpreal get_t();
  vector<string> get_data_string();
  
  mpreal get_energy();
  vector<mpreal> get_Ener();
  vector<mpreal> get_acc(int i);
};

#endif


