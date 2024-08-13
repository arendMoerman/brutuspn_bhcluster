#include "Bulirsch_Stoer.h"

Bulirsch_Stoer::Bulirsch_Stoer() {
  tolerance = "1e-6";
  n_max = 56;
  k_max = 128;
}
Bulirsch_Stoer::Bulirsch_Stoer(mpreal tolerance) {
  this->tolerance = tolerance;
  n_max = 56;
  k_max = 128;
}
Bulirsch_Stoer::Bulirsch_Stoer(mpreal tolerance, int n_max, int k_max) {
  this->tolerance = tolerance;
  this->n_max = n_max;
  this->k_max = k_max;
}

void Bulirsch_Stoer::set_tolerance(mpreal tolerance) {
  this->tolerance = tolerance;
}
void Bulirsch_Stoer::set_n_max(int n_max) {
  this->n_max = n_max;
}
void Bulirsch_Stoer::set_k_max(int k_max) {
  this->k_max = k_max;
}

mpreal Bulirsch_Stoer::get_tolerance() {
  return tolerance;
}
int Bulirsch_Stoer::get_n_max() {
  return n_max;
}
int Bulirsch_Stoer::get_k_max() {
  return k_max;
}

bool Bulirsch_Stoer::integrate(Cluster &cl, mpreal &dt) {
  Cluster c0 = cl;
  Cluster c  = cl;	
  mpreal timestep = dt;

  bool flag = step(c, timestep);
  if(flag) {
    cl = c;
    dt = timestep;
  }
  else {
    int k = 2;
    while( !flag && k <= k_max ) {
      timestep = dt/k;
      c = c0;
      flag = step(c, timestep);
      cout << "k = " << k << endl; 
      k += 2;
    }
  }  

  if(flag) {
    cl = c;
    dt = timestep;
  }

  return flag;
}
bool Bulirsch_Stoer::integrate(Cluster &cl, mpreal &t0, mpreal &t1) {
  int numStar = cl.s.size();
  x_sample.resize(numStar);
  y_sample.resize(numStar);
  z_sample.resize(numStar);
  vx_sample.resize(numStar);
  vy_sample.resize(numStar);
  vz_sample.resize(numStar);

  Cluster c0 = cl;
  Cluster c  = cl;	
  mpreal dt0 = t1-t0;
  mpreal timestep = dt0;

  bool flag = step(c, timestep);

  if(flag) {
    cl = c;
  }
  else {
    int k = 2;
    while( !flag && k <= k_max ) {
      timestep = dt0/k;
      c = c0;
      for(int q=0; q<k; q++) {
        flag = step(c, timestep);
        cout << "k = " << k << endl; 
        if(flag == false) break;
      }
      k += 2;
    }
  }  

  if(flag) {
    cl = c;
  }
  else {
    cerr << "Could not converge!" << endl;
    exit(1);
  }

  return flag;
}

bool Bulirsch_Stoer::step(Cluster &cl, mpreal &dt) {
  for (int i=0; i<x_sample.size(); i++) {
    x_sample[i].clear();
    y_sample[i].clear();
    z_sample[i].clear();
    vx_sample[i].clear();
    vy_sample[i].clear();
    vz_sample[i].clear();
  }

  bool flag;
  int n;
  int i = 0;
  vector<mpreal> h;
  vector<Cluster> c;
  Cluster cl_exp0 = cl;
  Cluster cl_exp  = cl;

  // n=1
  n=1;
  h.push_back( dt/n );
  c.push_back( cl );
  c[0].setw0();
  c[0].step( h[0] );

  for (int i = 0; i < x_sample.size(); i++) {
    x_sample[i].push_back(c[0].s[i].r[0]);
    y_sample[i].push_back(c[0].s[i].r[1]);
    z_sample[i].push_back(c[0].s[i].r[2]);
    vx_sample[i].push_back(c[0].s[i].v[0]);
    vy_sample[i].push_back(c[0].s[i].v[1]);
    vz_sample[i].push_back(c[0].s[i].v[2]);
  }

  cl_exp0 = c[0];

  // n=2
  n=2;
  h.push_back( dt/n );
  c.push_back( cl );
  c[1].setw0();
  c[1].step( h[1] );
  c[1].step( h[1] );

  extrapol( cl_exp, h, c );

  flag = error_control(cl_exp0, cl_exp);

  if(flag == false) {
    while(flag == false && n <= n_max) {

      n += 2;

      h.push_back( dt/n );
      c.push_back( cl );
      c[h.size()-1].setw0();
      for(int i=0; i<n; i++) c[h.size()-1].step( h[h.size()-1] ); //c[n/2].step( h[n/2] );

      cl_exp0 = cl_exp;
      extrapol( cl_exp, h, c );   

      flag = error_control(cl_exp0, cl_exp);
    }    
  }

  cl = cl_exp;
//cerr << n << " " << h.size() << endl;
  return flag;
}
void Bulirsch_Stoer::extrapol(Cluster &cl_exp, vector<mpreal> &dt, vector<Cluster> &c) {
  int N = c[0].s.size();
  int M = dt.size();

  for(int i=0; i<N; i++) {
    x_sample[i].push_back(  c[M-1].s[i].r[0] );
    y_sample[i].push_back(  c[M-1].s[i].r[1] );
    z_sample[i].push_back(  c[M-1].s[i].r[2] );
    vx_sample[i].push_back( c[M-1].s[i].v[0] );
    vy_sample[i].push_back( c[M-1].s[i].v[1] );
    vz_sample[i].push_back( c[M-1].s[i].v[2] );

    cl_exp.s[i].r[0] = extrapolate(dt, x_sample[i]);
    cl_exp.s[i].r[1] = extrapolate(dt, y_sample[i]);
    cl_exp.s[i].r[2] = extrapolate(dt, z_sample[i]);
    cl_exp.s[i].v[0] = extrapolate(dt, vx_sample[i]);
    cl_exp.s[i].v[1] = extrapolate(dt, vy_sample[i]);
    cl_exp.s[i].v[2] = extrapolate(dt, vz_sample[i]);
  }
}

// Neville's polynomial extrapolation algorithm
mpreal Bulirsch_Stoer::extrapolate(vector<mpreal> &x, vector<mpreal> &y) {
  int N = x.size();
  int i = N-1;
  for (int j=i-1; j >= 0; j--) {
    //y[j] = ( x[j]*y[j+1] -x[i]*y[j]  ) / ( x[j]-x[i] );

    mpreal f = x[j]/x[i];
    y[j] = y[j+1] + (y[j+1]-y[j])/(f*f-"1");
  }
  return y[0];
}
bool Bulirsch_Stoer::error_control(Cluster &c1, Cluster &c2) {
  int N = c1.s.size();
  bool flag = true;
  for(int i=0; i<N; i++) {
    if( fabs(c2.s[i].r[0]-c1.s[i].r[0]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].r[1]-c1.s[i].r[1]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].r[2]-c1.s[i].r[2]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[0]-c1.s[i].v[0]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[1]-c1.s[i].v[1]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[2]-c1.s[i].v[2]) > tolerance ) {
      flag = false;
      break;
    }
  }

  return flag;
}


