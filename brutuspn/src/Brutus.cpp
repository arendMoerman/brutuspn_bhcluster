#include "Brutus.h"

Brutus::Brutus(vector<mpreal> &data) {
  t = "0";
  this->data = data;
  N = data.size()/7;  

  tolerance = "1e-6";
  numBits = 56;

  eta = "0.2";

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  numBits = 56;

  eta = "0.2";

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  eta = "0.2";

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  this->eta = eta;

  setup();
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, int &nmax) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  this->eta = eta;

  setup(nmax);
}

mpreal Brutus::get_eta() {
  return eta;
}

void Brutus::setup() {
  Cluster c(data);
  cl = c;

  Bulirsch_Stoer b(tolerance);
  bs = b;
}
void Brutus::setup(int nmax) {
  Cluster c(data);
  cl = c;

  Bulirsch_Stoer b(tolerance, nmax, 128);
  bs = b;
}

bool Brutus::evolve(mpreal t_end) {
  bool merge;
  while (t<t_end) {
    cl.calcAcceleration_dt();

    dt = this->eta*cl.dt;

    mpreal t0 = t;
    t += dt;
    if(t > t_end) t = t_end;
    
    merge = cl.collisionDetection();
    if(merge) {
        N--;
        this->data = cl.get_data();
        return merge;
    }

    bool converged = bs.integrate(cl, t0, t);

    if(!converged) {
      cerr << "Not converged at " << t << "!" << endl;
      exit(1);
    }
  }
  this->data = cl.get_data();
  return merge;
}

void Brutus::reverse_velocities() {
    int N = cl.s.size();
    for(int i=0; i<N; i++) {
        for(int k=0; k<3; k++) {
            cl.s[i].v[k] *= "-1";
        }
    }
}

mpreal Brutus::get_t() {
  return t;
}
vector<mpreal> Brutus::get_data() {
  return data;
}
vector<double> Brutus::get_data_double() {
  int N = data.size()/7;
  vector<double> v(7*N, 0);
  for(int i=0; i<N; i++) {
    for(int j=0; j<7; j++) {
      v[i*7+j] = data[i*7+j].toDouble();
    }
  }
  return v;
}
vector<string> Brutus::get_data_string() {
  int N = data.size()/7;
  vector<string> v(7*N, "0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<7; j++) {
      v[i*7+j] = data[i*7+j].toString();
    }
  }
  return v;
}

mpreal Brutus::get_energy() {
    mpreal energy = cl.energies();
    return energy;
}

vector<mpreal> Brutus::get_Ener() {
    vector<mpreal> Ener = cl.get_Ener();
    return Ener;
}

vector<mpreal> Brutus::get_acc(int i) {
    int N = cl.s.size();
    
    vector<array<mpreal, 3>> ai;
    vector<mpreal> ao;
    ai = cl.getAcceleration();
    for(int k=0; k<N; k++) {
        ao.push_back(ai[i][k]);
    }
    return ao;
} 



