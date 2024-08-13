#include "Diagnostics.h"

// Conserved quantities

double Diagnostics::get_mass(vector<double> &data) {
  int N = data.size()/7;
  double M = 0;
  for(int i=0; i<N; i++) {
    M += data[i*7];
  }
  return M;
}
vector<double> Diagnostics::get_rcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> rcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rcm[j] += data[i*7]*data[i*7+(j+1)];
    }
  }
  for(int i=0; i<3; i++) rcm[i] /= M;
  return rcm;
}
vector<double> Diagnostics::get_vcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> vcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      vcm[j] += data[i*7]*data[i*7+(j+4)];
    }
  }
  for(int i=0; i<3; i++) vcm[i] /= M;
  return vcm;
}
vector<double> Diagnostics::get_lcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> lcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    lcm[0] += data[i*7]*(data[i*7+2]*data[i*7+6]-data[i*7+3]*data[i*7+5]);
    lcm[1] += data[i*7]*(data[i*7+3]*data[i*7+4]-data[i*7+1]*data[i*7+6]);
    lcm[2] += data[i*7]*(data[i*7+1]*data[i*7+5]-data[i*7+2]*data[i*7+4]);
  }
  for(int i=0; i<3; i++) lcm[i] /= M;
  return lcm;
}
double Diagnostics::get_kinetic_energy(vector<double> &data) {
  int N = data.size()/7;
  double ek = 0;
  for(int i=0; i<N; i++) {
    double m  = data[i*7];
    double vx = data[i*7+4];
    double vy = data[i*7+5];
    double vz = data[i*7+6];
    double v2 = vx*vx + vy*vy + vz*vz;
    ek += 0.5*m*v2;
  }
  return ek;
}
double Diagnostics::get_potential_energy(vector<double> &data) {
  int N = data.size()/7;
  double ep = 0;
  for(int i=0; i<N-1; i++) {
    double mi = data[i*7];
    double xi = data[i*7+1];
    double yi = data[i*7+2];
    double zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      double mj = data[j*7];
      double xj = data[j*7+1];
      double yj = data[j*7+2];
      double zj = data[j*7+3];

      double dx = xj - xi;
      double dy = yj - yi;
      double dz = zj - zi;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= mi*mj/sqrt(dr2);
    }
  }
  return ep;
}
double Diagnostics::get_energy(vector<double> &data) {
  double ek = get_kinetic_energy(data);
  double ep = get_potential_energy(data);
  return ek+ep;
}

// System properties

double Diagnostics::get_virial_radius(vector<double> &data) {
  double M = get_mass(data);
  double ep = get_potential_energy(data);
  double rv = -1*M*M / (2*ep);
  return rv;
}
double Diagnostics::get_harmonic_radius(vector<double> &data) {
  int N = data.size()/7;

  double ep = 0;
  for(int i=0; i<N-1; i++) {
    double xi = data[i*7+1];
    double yi = data[i*7+2];
    double zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      double xj = data[j*7+1];
      double yj = data[j*7+2];
      double zj = data[j*7+3];

      double dx = xj - xi;
      double dy = yj - yi;
      double dz = zj - zi;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= 1.0/sqrt(dr2);
    }
  }
  ep /= (N*N);

  double M = get_mass(data);
  double rh = -1*M*M / (2*ep);
  return rh;
}
double Diagnostics::get_velocity_disperion(vector<double> &data) {
  double M = get_mass(data);
  double ek = get_kinetic_energy(data);
  double sigma = sqrt(2*ek / M);
  return sigma;  
}

// Core collapse properties
vector<double> Diagnostics::get_rho(vector<double> &data) {
  int N = data.size()/7;

  vector<double> rho(N, 0);
  double pi = acos(-1.0);

  int numNN = 6;
  if(N < 8) numNN = N-1;

  Tools tools;
  vector< vector<int> > index_nn = tools.get_nn(data, numNN);
  
  for(int i=0; i<N; i++) {
    double mass = 0;
    for(int k=1; k<numNN-1; k++) {
      mass += data[index_nn[i][k]*7];
    }
    double dx = data[index_nn[i][numNN-1]*7+1] - data[i*7+1];
    double dy = data[index_nn[i][numNN-1]*7+2] - data[i*7+2];
    double dz = data[index_nn[i][numNN-1]*7+3] - data[i*7+3];
    double radius = sqrt(dx*dx + dy*dy + dz*dz);
    double volume = 4*pi/3 * radius*radius*radius;

    rho[i] = mass / volume;
  }

  return rho;
}
vector<double> Diagnostics::get_rho_center(vector<double> &data, vector<double> &rho) {
  vector<double> rho_center(3, 0);

  double rho_sum = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rho_center[j] += rho[i] * data[i*7+j];
    }
    rho_sum += rho[i];
  }  
  for(int j=0; j<3; j++) {
    rho_center[j] /= rho_sum;
  }

  return rho_center;
}
double Diagnostics::get_core_radius85(vector<double> &data, vector<double> &rho, vector<double> &rho_center) {
  double core_radius = 0;

  double rho_sum = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    double dx = data[i*7+1] - rho_center[0];
    double dy = data[i*7+2] - rho_center[1];
    double dz = data[i*7+3] - rho_center[2];
    double dr2 = dx*dx + dy*dy + dz*dz;
    core_radius += rho[i] * sqrt(dr2);
    rho_sum += rho[i];
  }
  core_radius /= rho_sum;

  return core_radius;
}
double Diagnostics::get_core_radius90(vector<double> &data, vector<double> &rho, vector<double> &rho_center) {
  double core_radius = 0;

  double rho_sum = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    double dx = data[i*7+1] - rho_center[0];
    double dy = data[i*7+2] - rho_center[1];
    double dz = data[i*7+3] - rho_center[2];
    double dr2 = dx*dx + dy*dy + dz*dz;
    core_radius += rho[i]*rho[i] * dr2;
    rho_sum += rho[i]*rho[i];
  }
  core_radius /= rho_sum;
  core_radius = sqrt(core_radius);

  return core_radius;
}
double Diagnostics::get_core_density(vector<double> &data, vector<double> &rho) {
  double core_density = 0;

  double rho_sum = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    core_density += rho[i] * rho[i];
    rho_sum += rho[i];
  }
  core_density /= rho_sum;

  return core_density;
}
double Diagnostics::get_core_mass(vector<double> &data, double &core_radius) {
  double core_mass = 0;

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    double r2 = 0;
    for(int j=0; j<3; j++) {
      r2 += data[i*7+j+1]*data[i*7+j+1];
    }
    if(r2 < core_radius*core_radius) {
      core_mass += data[i*7];
    }
  }

  return core_mass;
}
int Diagnostics::get_core_number(vector<double> &data, double &core_radius) {
  int core_number = 0;

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    double r2 = 0;
    for(int j=0; j<3; j++) {
      r2 += data[i*7+j+1]*data[i*7+j+1];
    }
    if(r2 < core_radius*core_radius) {
      core_number++;
    }
  }

  return core_number;
}






