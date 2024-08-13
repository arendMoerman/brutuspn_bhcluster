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

  Brutus(vector<mpreal> &data);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta);
  Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, int &nmax);

  mpreal get_eta();

  void setup();
  void setup(int nmax);

  bool evolve(mpreal t_end);
  void reverse_velocities();

  mpreal get_t();
  vector<mpreal> get_data();
  vector<double> get_data_double();
  vector<string> get_data_string();
  
  mpreal get_energy();
  vector<mpreal> get_Ener();
  vector<mpreal> get_acc(int i);
};

#endif


