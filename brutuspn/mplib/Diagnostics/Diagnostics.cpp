#include "Diagnostics.h"

// Conserved quantities

mpreal Diagnostics::get_mass(vector<mpreal> &data) {
  int N = data.size()/7;
  mpreal M = "0";
  for(int i=0; i<N; i++) {
    M += data[i*7];
  }
  return M;
}
vector<mpreal> Diagnostics::get_rcm(vector<mpreal> &data) {
  mpreal M = get_mass(data);
  vector<mpreal> rcm(3,"0");
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rcm[j] += data[i*7]*data[i*7+(j+1)];
    }
  }
  for(int i=0; i<3; i++) rcm[i] /= M;
  return rcm;
}
vector<mpreal> Diagnostics::get_vcm(vector<mpreal> &data) {
  mpreal M = get_mass(data);
  vector<mpreal> vcm(3,"0");
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      vcm[j] += data[i*7]*data[i*7+(j+4)];
    }
  }
  for(int i=0; i<3; i++) vcm[i] /= M;
  return vcm;
}
vector<mpreal> Diagnostics::get_lcm(vector<mpreal> &data) {
  mpreal M = get_mass(data);
  vector<mpreal> lcm(3,"0");
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    lcm[0] += data[i*7]*(data[i*7+2]*data[i*7+6]-data[i*7+3]*data[i*7+5]);
    lcm[1] += data[i*7]*(data[i*7+3]*data[i*7+4]-data[i*7+1]*data[i*7+6]);
    lcm[2] += data[i*7]*(data[i*7+1]*data[i*7+5]-data[i*7+2]*data[i*7+4]);
  }
  for(int i=0; i<3; i++) lcm[i] /= M;
  return lcm;
}
mpreal Diagnostics::get_kinetic_energy(vector<mpreal> &data) {
  int N = data.size()/7;
  mpreal ek = "0";
  for(int i=0; i<N; i++) {
    mpreal m  = data[i*7];
    mpreal vx = data[i*7+4];
    mpreal vy = data[i*7+5];
    mpreal vz = data[i*7+6];
    mpreal v2 = vx*vx + vy*vy + vz*vz;
    ek += "0.5"*m*v2;
  }
  return ek;
}
mpreal Diagnostics::get_potential_energy(vector<mpreal> &data) {
  int N = data.size()/7;
  mpreal ep = "0";
  for(int i=0; i<N-1; i++) {
    mpreal mi = data[i*7];
    mpreal xi = data[i*7+1];
    mpreal yi = data[i*7+2];
    mpreal zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      mpreal mj = data[j*7];
      mpreal xj = data[j*7+1];
      mpreal yj = data[j*7+2];
      mpreal zj = data[j*7+3];

      mpreal dx = xj - xi;
      mpreal dy = yj - yi;
      mpreal dz = zj - zi;
      mpreal dr2 = dx*dx + dy*dy + dz*dz;
      ep -= mi*mj/sqrt(dr2);
    }
  }
  return ep;
}
mpreal Diagnostics::get_energy(vector<mpreal> &data) {
  mpreal ek = get_kinetic_energy(data);
  mpreal ep = get_potential_energy(data);
  return ek+ep;
}

// System properties

mpreal Diagnostics::get_virial_radius(vector<mpreal> &data) {
  mpreal M = get_mass(data);
  mpreal ep = get_potential_energy(data);
  mpreal rv = "-1"*M*M / ("2"*ep);
  return rv;
}
mpreal Diagnostics::get_harmonic_radius(vector<mpreal> &data) {
  int N = data.size()/7;

  mpreal ep = "0";
  for(int i=0; i<N-1; i++) {
    mpreal xi = data[i*7+1];
    mpreal yi = data[i*7+2];
    mpreal zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      mpreal xj = data[j*7+1];
      mpreal yj = data[j*7+2];
      mpreal zj = data[j*7+3];

      mpreal dx = xj - xi;
      mpreal dy = yj - yi;
      mpreal dz = zj - zi;
      mpreal dr2 = dx*dx + dy*dy + dz*dz;
      ep -= "1"/sqrt(dr2);
    }
  }
  ep /= (N*N);

  mpreal M = get_mass(data);
  mpreal rh = "-1"*M*M / ("2"*ep);
  return rh;
}
mpreal Diagnostics::get_velocity_dispersion(vector<mpreal> &data) {
  int N = data.size()/7;
  mpreal Mtot = get_mass(data);
  vector<mpreal> vcm = get_vcm(data);

  mpreal ek = "0";
  for(int i=0; i<N; i++) {
    mpreal dvx = (data[i*7+4]-vcm[0]);
    mpreal dvy = (data[i*7+5]-vcm[1]);
    mpreal dvz = (data[i*7+6]-vcm[2]);
    mpreal dv2 = dvx*dvx + dvy*dvy + dvz*dvz;
    ek += "0.5"*data[i*7]*dv2;
  }

  mpreal sigma2 = "2"*ek/Mtot;  
  mpreal sigma = sqrt(sigma2);

  return sigma;  
}

vector<mpreal> Diagnostics::get_velocity_dispersion_1d(vector<mpreal> &data) {
  vector<mpreal> sigma(3, "0");

  int N = data.size()/7;
  mpreal Mtot = get_mass(data);
  vector<mpreal> vcm = get_vcm(data);
  mpreal ekx = "0";
  mpreal eky = "0";
  mpreal ekz = "0";
  for(int i=0; i<N; i++) {
    mpreal dvx = (data[i*7+4]-vcm[0]);
    mpreal dvy = (data[i*7+5]-vcm[1]);
    mpreal dvz = (data[i*7+6]-vcm[2]);

    ekx += "0.5"*data[i*7]*dvx*dvx;
    eky += "0.5"*data[i*7]*dvy*dvy;
    ekz += "0.5"*data[i*7]*dvz*dvz;
  }

  mpreal sigma2x = "2"*ekx/Mtot;  
  mpreal sigma2y = "2"*eky/Mtot;  
  mpreal sigma2z = "2"*ekz/Mtot;  

  sigma[0] = sqrt(sigma2x);
  sigma[1] = sqrt(sigma2y);
  sigma[2] = sqrt(sigma2z);

  return sigma;  
}

// Individual quantities

mpreal Diagnostics::get_kinetic_energy(vector<mpreal> &data, int index) {
  mpreal m  = data[index*7];
  mpreal vx = data[index*7+4];
  mpreal vy = data[index*7+5];
  mpreal vz = data[index*7+6];
  mpreal v2 = vx*vx + vy*vy + vz*vz;
  mpreal ek = "0.5"*m*v2;
  return ek;
}
mpreal Diagnostics::get_potential_energy(vector<mpreal> &data, int index) {
  int N = data.size()/7;
  mpreal ep = "0";

  mpreal mi = data[index*7];
  mpreal xi = data[index*7+1];
  mpreal yi = data[index*7+2];
  mpreal zi = data[index*7+3];
  for(int j=0; j<N; j++) {
    if(index != j) {
      mpreal mj = data[j*7];
      mpreal xj = data[j*7+1];
      mpreal yj = data[j*7+2];
      mpreal zj = data[j*7+3];

      mpreal dx = xj - xi;
      mpreal dy = yj - yi;
      mpreal dz = zj - zi;
      mpreal dr2 = dx*dx + dy*dy + dz*dz;
      ep -= mi*mj/sqrt(dr2);
    }
  }

  return ep;
}
mpreal Diagnostics::get_energy(vector<mpreal> &data, int index) {
  mpreal ek = get_kinetic_energy(data, index);
  mpreal ep = get_potential_energy(data, index);
  return ek+ep;
}

mpreal Diagnostics::get_radius(vector<mpreal> &data, int index) {
  mpreal x = data[index*7+1];
  mpreal y = data[index*7+2];
  mpreal z = data[index*7+3];
  mpreal R = sqrt(x*x + y*y + z*z);
  return R;
}
mpreal Diagnostics::get_speed(vector<mpreal> &data, int index) {
  mpreal vx = data[index*7+4];
  mpreal vy = data[index*7+5];
  mpreal vz = data[index*7+6];
  mpreal V = sqrt(vx*vx + vy*vy + vz*vz);
  return V;
}

// Core collapse properties

vector<mpreal> Diagnostics::get_rho(vector<mpreal> &data) {
  int N = data.size()/7;

  vector<mpreal> rho(N, "0");
  mpreal pi = acos("-1.0");

  int numNN = 6;
  if(N < 8) numNN = N-1;

  vector<double> ddata;
  for(int i=0; i<data.size(); i++) ddata.push_back(data[i].toDouble()); 

  Tools tools;
  vector< vector<int> > index_nn = tools.get_nn(ddata, numNN);
  
  for(int i=0; i<N; i++) {
    mpreal mass = "0";
    for(int k=1; k<numNN-1; k++) {
      mass += data[index_nn[i][k]*7];
    }
    mpreal dx = data[index_nn[i][numNN-1]*7+1] - data[i*7+1];
    mpreal dy = data[index_nn[i][numNN-1]*7+2] - data[i*7+2];
    mpreal dz = data[index_nn[i][numNN-1]*7+3] - data[i*7+3];
    mpreal radius = sqrt(dx*dx + dy*dy + dz*dz);
    mpreal volume = "4"*pi/"3" * radius*radius*radius;

    rho[i] = mass / volume;
  }

  return rho;
}
vector<mpreal> Diagnostics::get_rho_center(vector<mpreal> &data, vector<mpreal> &rho) {
  vector<mpreal> rho_center(3, "0");

  mpreal rho_sum = "0";
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rho_center[j] += rho[i] * data[i*7+j+1];
    }
    rho_sum += rho[i];
  }  
  for(int j=0; j<3; j++) {
    rho_center[j] /= rho_sum;
  }

  return rho_center;
}
mpreal Diagnostics::get_core_radius85(vector<mpreal> &data, vector<mpreal> &rho, vector<mpreal> &rho_center) {
  mpreal core_radius = "0";

  mpreal rho_sum = "0";
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1] - rho_center[0];
    mpreal dy = data[i*7+2] - rho_center[1];
    mpreal dz = data[i*7+3] - rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    core_radius += rho[i] * sqrt(dr2);
    rho_sum += rho[i];
  }
  core_radius /= rho_sum;

  return core_radius;
}
mpreal Diagnostics::get_core_radius90(vector<mpreal> &data, vector<mpreal> &rho, vector<mpreal> &rho_center) {
  mpreal core_radius = "0";

  mpreal rho_sum = "0";
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1] - rho_center[0];
    mpreal dy = data[i*7+2] - rho_center[1];
    mpreal dz = data[i*7+3] - rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    core_radius += rho[i]*rho[i] * dr2;
    rho_sum += rho[i]*rho[i];
  }
  core_radius /= rho_sum;
  core_radius = sqrt(core_radius);

  return core_radius;
}
mpreal Diagnostics::get_core_density(vector<mpreal> &data, vector<mpreal> &rho) {
  mpreal core_density = "0";

  mpreal rho_sum = "0";
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    core_density += rho[i] * rho[i];
    rho_sum += rho[i];
  }
  core_density /= rho_sum;

  return core_density;
}
mpreal Diagnostics::get_core_mass(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal &core_radius) {
  mpreal core_mass = "0";

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal r2 = "0";
    for(int j=0; j<3; j++) {
      mpreal x = data[i*7+j+1] - rho_center[j];
      r2 += x*x;
    }
    if(r2 < core_radius*core_radius) {
      core_mass += data[i*7];
    }
  }

  return core_mass;
}
int Diagnostics::get_core_number(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal &core_radius) {
  int core_number = 0;

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal r2 = "0";
    for(int j=0; j<3; j++) {
      mpreal x = data[i*7+j+1] - rho_center[j];
      r2 += x*x;
    }
    if(r2 < core_radius*core_radius) {
      core_number++;
    }
  }

  return core_number;
}

mpreal Diagnostics::get_halfmass_radius(vector<mpreal> &data, vector<mpreal> &rho_center) {
  int N = data.size()/7;

  vector<mpreal> R(N, "0");
  vector<mpreal> M(N, "0");
  mpreal Mtot = "0";

  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    
    R[i] = sqrt(dr2);
    M[i] = data[i*7];

    Mtot += data[i*7];
  }

  bool sorted = false;
  while(!sorted) {
    sorted = true;
    for(int i=0; i<N-1; i++) {
      if(R[i] > R[i+1]) {
        mpreal dummy = R[i];
        R[i] = R[i+1];
        R[i+1] = dummy;

        dummy = M[i];
        M[i] = M[i+1];
        M[i+1] = dummy;

        sorted = false;
      }
    }
  }

  mpreal Msum = "0";
  for(int i=0; i<N; i++) {
    Msum += M[i];
    if(Msum >= Mtot/"2") {
      return R[i];
    }
  }
}

mpreal Diagnostics::get_velocity_dispersion_halfmass(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rh) {
  int N = data.size()/7;

  mpreal Mtot = "0";
  vector<mpreal> vcm(3, "0");
  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    if(dr2 < Rh*Rh) {
      Mtot += data[i*7];
      vcm[0] += data[i*7]*data[i*7+4];
      vcm[1] += data[i*7]*data[i*7+5];
      vcm[2] += data[i*7]*data[i*7+6];
    }
  }
  for(int j=0; j<3; j++) {
    vcm[j] /= Mtot;
  }

  Mtot = "0";
  mpreal ek = "0";

  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    
    if(dr2 < Rh*Rh) {
      Mtot += data[i*7];

      mpreal dvx = data[i*7+4]-vcm[0];
      mpreal dvy = data[i*7+5]-vcm[1];
      mpreal dvz = data[i*7+6]-vcm[2];
      mpreal dv2 = dvx*dvx + dvy*dvy + dvz*dvz;        

      ek += "0.5"*data[i*7]*dv2;
    }
  }

  mpreal sigma2 = "2"*ek/Mtot;
  mpreal sigma = sqrt(sigma2);

  return sigma;
}

vector<mpreal> Diagnostics::get_velocity_dispersion_halfmass_1d(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rh) {
  vector<mpreal> sigma(3, "0");

  int N = data.size()/7;

  mpreal Mtot = "0";
  vector<mpreal> vcm(3, "0");
  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;

    if(dr2 < Rh*Rh) {
      Mtot += data[i*7];
      vcm[0] += data[i*7]*data[i*7+4];
      vcm[1] += data[i*7]*data[i*7+5];
      vcm[2] += data[i*7]*data[i*7+6];
    }
  }

  for(int j=0; j<3; j++) {
    vcm[j] /= Mtot;
  }

  Mtot = "0";
  mpreal ekx = "0";
  mpreal eky = "0";
  mpreal ekz = "0";

  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    
    if(dr2 < Rh*Rh) {
      Mtot += data[i*7];

      mpreal dvx = data[i*7+4]-vcm[0];
      mpreal dvy = data[i*7+5]-vcm[1];
      mpreal dvz = data[i*7+6]-vcm[2];

      ekx += "0.5"*data[i*7]*dvx*dvx;
      eky += "0.5"*data[i*7]*dvy*dvy;
      ekz += "0.5"*data[i*7]*dvz*dvz;
    }
  }

  mpreal sigma2x = "2"*ekx/Mtot;
  mpreal sigma2y = "2"*eky/Mtot;
  mpreal sigma2z = "2"*ekz/Mtot;

  sigma[0] = sqrt(sigma2x);
  sigma[1] = sqrt(sigma2y);
  sigma[2] = sqrt(sigma2z);

  return sigma;
}

vector<mpreal> Diagnostics::get_velocity_dispersion_plummer_radius_1d(vector<mpreal> &data, vector<mpreal> &rho_center, mpreal Rp) {
  vector<mpreal> sigma(3, "0");

  int N = data.size()/7;

  mpreal Mtot = "0";
  vector<mpreal> vcm(3, "0");
  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;

    if(dr2 < Rp*Rp) {
      Mtot += data[i*7];
      vcm[0] += data[i*7]*data[i*7+4];
      vcm[1] += data[i*7]*data[i*7+5];
      vcm[2] += data[i*7]*data[i*7+6];
    }
  }

  for(int j=0; j<3; j++) {
    vcm[j] /= Mtot;
  }

  Mtot = "0";
  mpreal ekx = "0";
  mpreal eky = "0";
  mpreal ekz = "0";

  for(int i=0; i<N; i++) {
    mpreal dx = data[i*7+1]-rho_center[0];
    mpreal dy = data[i*7+2]-rho_center[1];
    mpreal dz = data[i*7+3]-rho_center[2];
    mpreal dr2 = dx*dx + dy*dy + dz*dz;
    
    if(dr2 < Rp*Rp) {
      Mtot += data[i*7];

      mpreal dvx = data[i*7+4]-vcm[0];
      mpreal dvy = data[i*7+5]-vcm[1];
      mpreal dvz = data[i*7+6]-vcm[2];

      ekx += "0.5"*data[i*7]*dvx*dvx;
      eky += "0.5"*data[i*7]*dvy*dvy;
      ekz += "0.5"*data[i*7]*dvz*dvz;
    }
  }

  mpreal sigma2x = "2"*ekx/Mtot;
  mpreal sigma2y = "2"*eky/Mtot;
  mpreal sigma2z = "2"*ekz/Mtot;

  sigma[0] = sqrt(sigma2x);
  sigma[1] = sqrt(sigma2y);
  sigma[2] = sqrt(sigma2z);

  return sigma;
}

vector<mpreal> Diagnostics::get_radii(vector<mpreal> &data, vector<mpreal> &rho_center) {
  int N = data.size()/7;

  vector<mpreal> rs;

  for(int i=0; i<N; i++) {

    mpreal x = data[i*7+1]-rho_center[0];
    mpreal y = data[i*7+2]-rho_center[1];
    mpreal z = data[i*7+3]-rho_center[2];
    mpreal R = sqrt(x*x + y*y + z*z);

    rs.push_back(R);
  }

  return rs;
}





