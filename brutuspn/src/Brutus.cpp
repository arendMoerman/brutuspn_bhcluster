#include "Brutus.h"

Brutus::Brutus(vector<mpreal> &data, mpreal r_merge) {
  t = "0";
  this->data = data;
  N = data.size()/7;  

  tolerance = "1e-6";
  numBits = 56;

  eta = "0.2";

  setup(r_merge);
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, mpreal r_merge) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  numBits = 56;

  eta = "0.2";

  setup(r_merge);
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal r_merge) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  eta = "0.2";

  setup(r_merge);
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, mpreal r_merge) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  this->eta = eta;

  setup(r_merge);
}
Brutus::Brutus(mpreal &t, vector<mpreal> &data, mpreal &tolerance, int &numBits, mpreal &eta, int &nmax, mpreal r_merge) {
  this->t = t;
  this->data = data;
  N = data.size()/7;  

  this->tolerance = tolerance;
  this->numBits = numBits;

  this->eta = eta;

  setup(nmax, r_merge);
}

mpreal Brutus::get_eta() {
  return eta;
}

void Brutus::setup(mpreal r_merge) {
  Cluster c(data, r_merge);
  cl = c;

  Bulirsch_Stoer b(tolerance);
  bs = b;
}
void Brutus::setup(int nmax, mpreal r_merge) {
  Cluster c(data, r_merge);
  cl = c;

  Bulirsch_Stoer b(tolerance, nmax, 128);
  bs = b;
}

void Brutus::evolve(mpreal t_end, mergerOut &merge) {
  while (t<t_end) {
    cl.calcAcceleration_dt();

    dt = this->eta*cl.dt;

    mpreal t0 = t;
    t += dt;
    if(t > t_end) t = t_end;
    
    cl.collisionDetection(merge);
    if(merge.merged) {
        this->data = cl.get_data();
        break;
    }

    bool converged = bs.integrate(cl, t0, t);

    if(!converged) {
      cerr << "Not converged at " << t << "!" << endl;
      exit(1);
    }
  }
  this->data = cl.get_data();
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

vector<string> Brutus::get_data_string() {
  int N = data.size()/7;
  vector<string> v(13*N, "0");
  
  vector<array<mpreal, 3>> acc;
  vector<array<mpreal, 3>> jerk;

  acc = cl.getAcceleration();   
  jerk = cl.getJerk();   

  for(int i=0; i<N; i++) {
    for(int j=0; j<7; j++) {
      v[i*13+j] = data[i*7+j].toString();
    }
    
    for(int j=0; j<3; j++) {
      v[i*13+7+j] = acc[i][j].toString();
      v[i*13+10+j] = jerk[i][j].toString();
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



