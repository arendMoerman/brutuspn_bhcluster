#include "Products.h"

/////////////////////////////////////////////////////////////////////

Products::Products() {
  Resc = 4;
  fiso_bin = 4;
  fiso_flyby = 4;
}
Products::Products(double fiso) {
  Resc = 4;
  fiso_bin = fiso;
  fiso_flyby = fiso;
}

/////////////////////////////////////////////////////////////////////

vector<double> Products::get_cm_multiple(vector<double> &data, vector<int> &index) {
  int M = index.size();

  vector<double> cm(7, 0);

  for(int i=0; i<M; i++) {
    cm[0] += data[index[i]*7+0];
    cm[1] += data[index[i]*7+0]*data[index[i]*7+1];
    cm[2] += data[index[i]*7+0]*data[index[i]*7+2];
    cm[3] += data[index[i]*7+0]*data[index[i]*7+3];
    cm[4] += data[index[i]*7+0]*data[index[i]*7+4];
    cm[5] += data[index[i]*7+0]*data[index[i]*7+5];
    cm[6] += data[index[i]*7+0]*data[index[i]*7+6];
  }
  for(int i=1; i<7; i++) {
    cm[i] /= cm[0];
  }

  return cm;
}
vector<double> Products::get_coordinates_multiple(vector<double> &data, vector<int> &index) {
  vector<double> multiple_data;

  int M = index.size();
  for(int i=0; i<M; i++) {
    int myindex = index[i];
    for(int j=0; j<7; j++) {
      multiple_data.push_back(data[myindex*7+j]);
    }
  }

  return multiple_data;
}
vector<double> Products::get_coordinates_multiple_in_cm_frame(vector<double> &data, vector<int> &index) {
  vector<double> multiple_data = get_coordinates_multiple(data, index);
  vector<double> cm = get_cm_multiple(data, index);

  int M = index.size();
  for(int i=0; i<M; i++) {
    for(int j=1; j<7; j++) {
      multiple_data[i*7+j] -= cm[0*7+j];
    }
  }

  return multiple_data;
}

double Products::get_kinetic_energy_particle(vector<double> &data, int index) {
  double m  = data[index*7];
  double vx = data[index*7+4];
  double vy = data[index*7+5];
  double vz = data[index*7+6];
  double v2 = vx*vx + vy*vy + vz*vz;
  double ek = 0.5*m*v2;
  return ek;
}
double Products::get_potential_energy_particle(vector<double> &data, int index) {
  double m  = data[index*7];
  double x  = data[index*7+1];
  double y  = data[index*7+2];
  double z  = data[index*7+3];
  double ep = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    if(i != index) {
      double mi  = data[i*7];
      double xi  = data[i*7+1];
      double yi  = data[i*7+2];
      double zi  = data[i*7+3];
      double dx = xi - x;
      double dy = yi - y;
      double dz = zi - z;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= m*mi / sqrt(dr2);
    }
  }  
  return ep;
}
double Products::get_energy_particle(vector<double> &data, int index) {
  double ek = get_kinetic_energy_particle(data, index);
  double ep = get_potential_energy_particle(data, index);
  return ek+ep;
}

double Products::get_external_kinetic_energy_multiple(vector<double> &data, vector<int> index) {
  vector<double> cm = get_cm_multiple(data, index);
  double m  = cm[0];
  double vx = cm[4];
  double vy = cm[5];
  double vz = cm[6];
  double v2 = vx*vx + vy*vy + vz*vz;
  double ek = 0.5*m*v2;
  return ek;
}
double Products::get_external_potential_energy_multiple(vector<double> &data, vector<int> index) {
  int M = index.size();

  vector<double> cm = get_cm_multiple(data, index);
  double m  = cm[0];
  double x  = cm[1];
  double y  = cm[2];
  double z  = cm[3];

  double ep = 0;
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double mi  = data[i*7];
      double xi  = data[i*7+1];
      double yi  = data[i*7+2];
      double zi  = data[i*7+3];
      double dx = xi - x;
      double dy = yi - y;
      double dz = zi - z;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= m*mi / sqrt(dr2);
    }
  }  
  return ep;
}
double Products::get_external_energy_multiple(vector<double> &data, vector<int> index) {
  double ek = get_external_kinetic_energy_multiple(data, index);
  double ep = get_external_potential_energy_multiple(data, index);
  return ek+ep;
}

double Products::get_internal_kinetic_energy_multiple(vector<double> &data, vector<int> index) {
  int M = index.size();

  vector<double> cm = get_cm_multiple(data, index);
  double vxcm = cm[4];
  double vycm = cm[5];
  double vzcm = cm[6];

  double ek = 0;
  for(int i=0; i<M; i++) {
    double m  = data[index[i]*7];
    double vx = data[index[i]*7+4]-vxcm;
    double vy = data[index[i]*7+5]-vycm;
    double vz = data[index[i]*7+6]-vzcm;
    double v2 = vx*vx+vy*vy+vz*vz;
    ek += 0.5*m*v2;
  }

  return ek;
}
double Products::get_internal_potential_energy_multiple(vector<double> &data, vector<int> index) {
  int M = index.size();

  double ep = 0;
  for(int i=0; i<M-1; i++) {
    double m1 = data[index[i]*7];
    double x1 = data[index[i]*7+1];
    double y1 = data[index[i]*7+2];
    double z1 = data[index[i]*7+3];
    for(int j=i+1; j<M; j++) {
      double m2 = data[index[j]*7];
      double x2 = data[index[j]*7+1];
      double y2 = data[index[j]*7+2];
      double z2 = data[index[j]*7+3];
      
      double dx = x2-x1;
      double dy = y2-y1;
      double dz = z2-z1;
      double dr2 = dx*dx+dy*dy+dz*dz;
      
      ep -= m1*m2/sqrt(dr2);
    }
  }

  return ep;
}
double Products::get_internal_energy_multiple(vector<double> &data, vector<int> index) {
  double ek = get_internal_kinetic_energy_multiple(data, index);
  double ep = get_internal_potential_energy_multiple(data, index);
  return ek+ep;
}

double Products::get_angular_momentum_particle(vector<double> &data, int index) {
  double m = data[index*7+0];
  double x = data[index*7+1];
  double y = data[index*7+2];
  double z = data[index*7+3];
  double vx = data[index*7+4];
  double vy = data[index*7+5];
  double vz = data[index*7+6];
  
  double lx = y*vz-z*vy;
  double ly = z*vx-x*vz;
  double lz = x*vy-y*vx;
  double l2 = lx*lx + ly*ly + lz*lz;

  double l = m*sqrt(l2);

  return l;
}
double Products::get_external_angular_momentum_multiple(vector<double> &data, vector<int> &index) { 
  vector<double> cm = get_cm_multiple(data, index);

  double m  = cm[0];
  double x  = cm[1];
  double y  = cm[2];
  double z  = cm[3];
  double vx = cm[4];
  double vy = cm[5];
  double vz = cm[6];

  double lx = y*vz-z*vy;
  double ly = z*vx-x*vz;
  double lz = x*vy-y*vx;
  double l2 = lx*lx + ly*ly + lz*lz;

  double l = m*sqrt(l2);

  return l;
}
double Products::get_internal_angular_momentum_multiple(vector<double> &data, vector<int> &index) {
  int M = index.size();

  vector<double> cm = get_cm_multiple(data, index);
  
  vector<double> subset(7*M, 0);
  for(int i=0; i<M; i++) {
    subset[i*7+0] = data[index[i]*7+0];
    for(int j=1; j<7; j++) {
      subset[i*7+j] = data[index[i]*7+j] - cm[j];
    }
  }

  double Lx=0, Ly=0, Lz=0;
  for(int i=0; i<M; i++) {
    double m  = subset[i*7+0];
    double x  = subset[i*7+1];
    double y  = subset[i*7+2];
    double z  = subset[i*7+3];
    double vx = subset[i*7+4];
    double vy = subset[i*7+5];
    double vz = subset[i*7+6];

    double lx = m*(y*vz-z*vy);
    double ly = m*(z*vx-x*vz);
    double lz = m*(x*vy-y*vx);

    Lx += lx;
    Ly += ly;
    Lz += lz;
  }

  double L2 = Lx*Lx + Ly*Ly + Lz*Lz;
  double L = sqrt(L2);

  return L;
}

/////////////////////////////////////////////////////////////////////

bool Products::isIsolated_abs_particle(vector<double> &data, int index, double R) {
  int N = data.size()/7;

  double xcm = data[index*7+1]; 
  double ycm = data[index*7+2]; 
  double zcm = data[index*7+3]; 

  double R2 = R*R;

  for(int i=0; i<N; i++) {
    if(i != index) {
      double dx = data[i*7+1]-xcm;
      double dy = data[i*7+2]-ycm;
      double dz = data[i*7+3]-zcm;
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 < R2) return false;
    }
  }

  return true;
}
bool Products::isIsolated_abs_multiple(vector<double> &data, vector<int> &index, double R) {
  int N = data.size()/7;
  int M = index.size();

  double mcm=0, xcm=0, ycm=0, zcm=0;
  for(int i=0; i<M; i++) {
    mcm += data[index[i]*7];
    xcm += data[index[i]*7]*data[index[i]*7+1];
    ycm += data[index[i]*7]*data[index[i]*7+2];
    zcm += data[index[i]*7]*data[index[i]*7+3];
  }
  xcm /= mcm;
  ycm /= mcm;
  zcm /= mcm;

  double R2 = R*R;

  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double dx = data[i*7+1]-xcm;
      double dy = data[i*7+2]-ycm;
      double dz = data[i*7+3]-zcm;
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 < R2) {
        return false;
      }	
    }
  }

  return true;
}
bool Products::isIsolated_rel_multiple(vector<double> &data, vector<int> &index, double f) {
  int N = data.size()/7;
  int M = index.size();

  double mcm=0, xcm=0, ycm=0, zcm=0;
  for(int i=0; i<M; i++) {
    mcm += data[index[i]*7];
    xcm += data[index[i]*7]*data[index[i]*7+1];
    ycm += data[index[i]*7]*data[index[i]*7+2];
    zcm += data[index[i]*7]*data[index[i]*7+3];
  }
  xcm /= mcm;
  ycm /= mcm;
  zcm /= mcm;

  double f2 = f*f;

  double R2 = 0;
  for(int i=0; i<M; i++) {
    double dx = data[index[i]*7+1]-xcm;
    double dy = data[index[i]*7+2]-ycm;
    double dz = data[index[i]*7+3]-zcm;
    double dr2 = dx*dx+dy*dy+dz+dz;
    if(dr2 > R2) R2 = dr2;
  }

  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double dx = data[i*7+1]-xcm;
      double dy = data[i*7+2]-ycm;
      double dz = data[i*7+3]-zcm;
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 < f2*R2) {
        return false;
      }	
    }
  }

  return true;
}

/////////////////////////////////////////////////////////////////////

bool Products::is_nn_multiple(vector<double> &data, vector<int> index) {
  int N = data.size()/7;
  int M = index.size();

  double R2 = 0;
  for(int i=0; i<M-1; i++) {
    for(int j=i+1; j<M; j++) {
      double dx = data[index[i]*7+1]-data[index[j]*7+1];
      double dy = data[index[i]*7+2]-data[index[j]*7+2];
      double dz = data[index[i]*7+3]-data[index[j]*7+3];
      double dr2 = dx*dx + dy*dy + dz*dz;
      if(dr2 > R2) R2 = dr2;
    }
  }

  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double myR2 = 1e100;
      for(int j=0; j<M; j++) {
        double dx = data[i*7+1]-data[index[j]*7+1];
        double dy = data[i*7+2]-data[index[j]*7+2];
        double dz = data[i*7+3]-data[index[j]*7+3];
        double dr2 = dx*dx + dy*dy + dz*dz;
        if(dr2 < myR2) myR2 = dr2;
      }
      if(myR2 < R2) {
        return false;
      }
    }
  }

  return true;
}

/////////////////////////////////////////////////////////////////////

vector<int> Products::get_nn_list_particle(vector<double> &data, int index, int numNN) {
  int N = data.size()/7;

  double x0 = data[index*7+1];
  double y0 = data[index*7+2];
  double z0 = data[index*7+3];

  vector<int> nn_index(numNN, 0);
  vector<double> nn_dr2(numNN, 1e100);

  for(int i=0; i<N; i++) {
    if(i != index) {
      double dx = data[i*7+1]-x0;
      double dy = data[i*7+2]-y0;
      double dz = data[i*7+3]-z0;
      double dr2 = dx*dx + dy*dy + dz*dz;  
      for(int j=0; j<numNN; j++) {
        if(dr2 < nn_dr2[j]) {
          for(int k=numNN-1; k>j; k--) {
	    nn_dr2[k] = nn_dr2[k-1];
	    nn_index[k] = nn_index[k-1];
          }
          nn_dr2[j] = dr2;
          nn_index[j] = i;
          break;
        }  
      }
    }
  }

  return nn_index;
}
vector<int> Products::get_nn_list_multiple(vector<double> &data, vector<int> index, int numNN) {
  int N = data.size()/7;
  int M = index.size();

  double m0=0, x0=0, y0=0, z0=0;
  for(int i=0; i<M; i++) {
    m0 += data[index[i]*7];
    x0 += data[index[i]*7]*data[index[i]*7+1];
    y0 += data[index[i]*7]*data[index[i]*7+2];
    z0 += data[index[i]*7]*data[index[i]*7+3];
  }
  x0 /= m0;
  y0 /= m0;
  z0 /= m0;

  vector<int> nn_index(numNN, 0);
  vector<double> nn_dr2(numNN, 1e100);

  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double dx = data[i*7+1]-x0;
      double dy = data[i*7+2]-y0;
      double dz = data[i*7+3]-z0;
      double dr2 = dx*dx + dy*dy + dz*dz;  
      for(int j=0; j<numNN; j++) {
        if(dr2 < nn_dr2[j]) {
          for(int k=numNN-1; k>j; k--) {
	    nn_dr2[k] = nn_dr2[k-1];
	    nn_index[k] = nn_index[k-1];
          }
          nn_dr2[j] = dr2;
          nn_index[j] = i;
          break;
        }  
      }
    }
  }

  return nn_index;
}

double Products::get_distance_nn_particle(vector<double> &data, int index) {
  int N = data.size()/7;

  double x0 = data[index*7+1];
  double y0 = data[index*7+2];
  double z0 = data[index*7+3];

  double dr2_min = 1e100;

  for(int i=0; i<N; i++) {
    if(i != index) {
      double dx = data[i*7+1]-x0;
      double dy = data[i*7+2]-y0;
      double dz = data[i*7+3]-z0;
      double dr2 = dx*dx + dy*dy + dz*dz;  
      if(dr2 < dr2_min) dr2_min = dr2;
    }
  }

  double dr_min = sqrt(dr2_min);

  return dr_min;
}
double Products::get_distance_nn_multiple(vector<double> &data, vector<int> index) {
  int N = data.size()/7;
  int M = index.size();

  double m0=0, x0=0, y0=0, z0=0;
  for(int i=0; i<M; i++) {
    m0 += data[index[i]*7];
    x0 += data[index[i]*7]*data[index[i]*7+1];
    y0 += data[index[i]*7]*data[index[i]*7+2];
    z0 += data[index[i]*7]*data[index[i]*7+3];
  }
  x0 /= m0;
  y0 /= m0;
  z0 /= m0;

  double dr2_min = 1e100;

  for(int i=0; i<N; i++) {
    bool original = true;
    for(int j=0; j<M; j++) {
      if(index[j] == i) {
	original = false;
	break;
      }  
    }
    if(original) {
      double dx = data[i*7+1]-x0;
      double dy = data[i*7+2]-y0;
      double dz = data[i*7+3]-z0;
      double dr2 = dx*dx + dy*dy + dz*dz;  
      if(dr2 < dr2_min) dr2_min = dr2;      
    }
  }

  double dr_min = sqrt(dr2_min);

  return dr_min;
}

/////////////////////////////////////////////////////////////////////

vector<double> Products::get_escaper_properties(vector<double> &data, int index1) {
// e, l2, m, r, v, rdotv, theta, phi

  double m  = data[index1*7];
  double x  = data[index1*7+1];
  double y  = data[index1*7+2];
  double z  = data[index1*7+3];
  double vx = data[index1*7+4];
  double vy = data[index1*7+5];
  double vz = data[index1*7+6];

  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double v2 = vx*vx + vy*vy + vz*vz;
  double v = sqrt(v2);

  double rdotv = x*vx + y*vy + z*vz;

  double theta = acos(z/r);
  double phi = atan2(y, x);

  double e = get_energy_particle(data, index1);
  double l = get_angular_momentum_particle(data, index1);

  double l2 = l*l;

  vector<double> orbelem(8);
  orbelem[0] = e; 
  orbelem[1] = l2;
  orbelem[2] = m;
  orbelem[3] = r;
  orbelem[4] = v;
  orbelem[5] = rdotv;
  orbelem[6] = theta;
  orbelem[7] = phi;

  return orbelem;
}

vector<double> Products::get_binary_properties(vector<double> &data, int index1, int index2) {
// a, e, EB, E_int, L2_int, mu, Rcm, Vcm, RcmdotVcm, E_ext, rp, theta, rdotv, r

  double m1  = data[index1*7];
  double x1  = data[index1*7+1];
  double y1  = data[index1*7+2];
  double z1  = data[index1*7+3];
  double vx1 = data[index1*7+4];
  double vy1 = data[index1*7+5];
  double vz1 = data[index1*7+6];

  double m2  = data[index2*7];
  double x2  = data[index2*7+1];
  double y2  = data[index2*7+2];
  double z2  = data[index2*7+3];
  double vx2 = data[index2*7+4];
  double vy2 = data[index2*7+5];
  double vz2 = data[index2*7+6];

  double mu = m1+m2;
  double xcm = (m1*x1 + m2*x2)/mu;
  double ycm = (m1*y1 + m2*y2)/mu;
  double zcm = (m1*z1 + m2*z2)/mu;
  double vxcm = (m1*vx1 + m2*vx2)/mu;
  double vycm = (m1*vy1 + m2*vy2)/mu;
  double vzcm = (m1*vz1 + m2*vz2)/mu;
  double r2cm = xcm*xcm + ycm*ycm + zcm*zcm;
  double rcm = sqrt(r2cm);
  double v2cm = vxcm*vxcm + vycm*vycm + vzcm*vzcm;  
  double vcm = sqrt(v2cm);
  double rcmdotvcm = xcm*vxcm + ycm*vycm + zcm*vzcm;

  vector<int> index;
  index.push_back(index1);
  index.push_back(index2);
  double E_ext = get_external_energy_multiple(data, index);

  double x = x2-x1;
  double y = y2-y1;
  double z = z2-z1;
  double vx = vx2-vx1;
  double vy = vy2-vy1;
  double vz = vz2-vz1;  
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double v2 = vx*vx + vy*vy + vz*vz;
  double rdotv = x*vx + y*vy + z*vz;

  double e = 0.5*v2 - mu/r;
  double lx = y*vz-z*vy;
  double ly = z*vx-x*vz;
  double lz = x*vy-y*vx;  
  double l2 = lx*lx+ly*ly+lz*lz;

  double mya = -1;
  double mye = -1;
  double myEb = 0;
  double rp = 0;
  double theta = 0; 

  if(e < 0) {
    double D = mu*mu + 2.0*e*l2;
    if(D >= 0) {
      rp = (-mu+sqrt(D)) / (2.0*e);  
      double ra = (-mu-sqrt(D)) / (2.0*e);    
      mya  = 0.5*(rp+ra);
      mye  = 1.0-2.0/(ra/rp+1.0);
      myEb = 0.5*m1*m2/mya;
      //theta = acos( ((mya*(1.0-mye*mye) / r - 1.0))/mye );
    }
  }
  else if(e == 0) {
    rp = l2/(2.0*mu);
    //mya   = ;
    //theta = acos( l2/mu/r-1.0 );
  }
  else {
    double D = mu*mu + 2.0*e*l2;
    if(D >= 0) {
      rp = (-mu+sqrt(D)) / (2.0*e);  
      //mya   = ;
      //mye   = ;
      //theta = acos( (l2/mu/r-1.0)/mye );
    }  
  }

  vector<double> orbelem(14);
  orbelem[0] = mya; 
  orbelem[1] = mye;
  orbelem[2] = myEb;
  orbelem[3] = e;
  orbelem[4] = l2;
  orbelem[5] = mu;
  orbelem[6] = rcm;
  orbelem[7] = vcm;
  orbelem[8] = rcmdotvcm;
  orbelem[9] = E_ext;
  orbelem[10] = rp;
  orbelem[11] = theta;
  orbelem[12] = rdotv;
  orbelem[13] = r;

  return orbelem;
}

vector<double> Products::get_triple_properties(vector<double> &data, int index1, int index2, int index3) {
// a1, e1, Eb1, a2, e2, Eb2, f_stab, theta_int, E_int, L2_int, mu, Rcm, Vcm, RcmdotVcm, E_ext

  double m1  = data[index1*7];
  double x1  = data[index1*7+1];
  double y1  = data[index1*7+2];
  double z1  = data[index1*7+3];
  double vx1 = data[index1*7+4];
  double vy1 = data[index1*7+5];
  double vz1 = data[index1*7+6];

  double m2  = data[index2*7];
  double x2  = data[index2*7+1];
  double y2  = data[index2*7+2];
  double z2  = data[index2*7+3];
  double vx2 = data[index2*7+4];
  double vy2 = data[index2*7+5];
  double vz2 = data[index2*7+6];

  double m3  = data[index3*7];
  double x3  = data[index3*7+1];
  double y3  = data[index3*7+2];
  double z3  = data[index3*7+3];
  double vx3 = data[index3*7+4];
  double vy3 = data[index3*7+5];
  double vz3 = data[index3*7+6];

  double mu = m1+m2+m3;
  double xcm = (m1*x1 + m2*x2 + m3*x3)/mu;
  double ycm = (m1*y1 + m2*y2 + m3*y3)/mu;
  double zcm = (m1*z1 + m2*z2 + m3*z3)/mu;
  double vxcm = (m1*vx1 + m2*vx2 + m3*vx3)/mu;
  double vycm = (m1*vy1 + m2*vy2 + m3*vy3)/mu;
  double vzcm = (m1*vz1 + m2*vz2 + m3*vz3)/mu;

  double r2cm = xcm*xcm + ycm*ycm + zcm*zcm;
  double rcm = sqrt(r2cm);
  double v2cm = vxcm*vxcm + vycm*vycm + vzcm*vzcm;  
  double vcm = sqrt(v2cm);

  double rcmdotvcm = xcm*vxcm + ycm*vycm + zcm*vzcm;

  vector<int> index;
  index.push_back(index1);
  index.push_back(index2);
  index.push_back(index3);
  double E_ext = get_external_energy_multiple(data, index);

  double e = get_internal_energy_multiple(data, index);
  double l = get_internal_angular_momentum_multiple(data, index);
  double l2 = l*l;

  double dr2_12 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
  double dr2_13 = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1);
  double dr2_23 = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2);
  int indexa = index1, indexb = index2, indexc = index3;
  if(dr2_13 < dr2_12 && dr2_13 < dr2_23) {
    indexa = index1;
    indexb = index3;
    indexc = index2; 
  }
  if(dr2_23 < dr2_12 && dr2_23 < dr2_13) {
    indexa = index2;
    indexb = index3;
    indexc = index1; 
  }

  vector<double> orbelem_inner = get_binary_properties(data, indexa, indexb); 

  double m_inner = data[indexa*7] + data[indexb*7];
  double x_inner = (data[indexa*7]*data[indexa*7+1] + data[indexb*7]*data[indexb*7+1])/m_inner;  
  double y_inner = (data[indexa*7]*data[indexa*7+2] + data[indexb*7]*data[indexb*7+2])/m_inner;
  double z_inner = (data[indexa*7]*data[indexa*7+3] + data[indexb*7]*data[indexb*7+3])/m_inner;
  double vx_inner = (data[indexa*7]*data[indexa*7+4] + data[indexb*7]*data[indexb*7+4])/m_inner;  
  double vy_inner = (data[indexa*7]*data[indexa*7+5] + data[indexb*7]*data[indexb*7+5])/m_inner;
  double vz_inner = (data[indexa*7]*data[indexa*7+6] + data[indexb*7]*data[indexb*7+6])/m_inner;

  vector<double> data2;
  data2.push_back(m_inner);
  data2.push_back(x_inner);
  data2.push_back(y_inner);
  data2.push_back(z_inner);
  data2.push_back(vx_inner);
  data2.push_back(vy_inner);
  data2.push_back(vz_inner);
  data2.push_back(data[indexc*7+0]);
  data2.push_back(data[indexc*7+1]);
  data2.push_back(data[indexc*7+2]);
  data2.push_back(data[indexc*7+3]);
  data2.push_back(data[indexc*7+4]);
  data2.push_back(data[indexc*7+5]);
  data2.push_back(data[indexc*7+6]);

  vector<double> orbelem_outer = get_binary_properties(data2, 0, 1); 

  double a1 = orbelem_inner[0];
  double e1 = orbelem_inner[1];
  double Eb1 = orbelem_inner[2];
  double a2 = orbelem_outer[0];
  double e2 = orbelem_outer[1];
  double Eb2 = orbelem_outer[2];

  double rp = a2*(1.0-e2);
  double f_stab = rp/a1;

  double dxi  = data[indexb*7+1]-data[indexa*7+1];
  double dyi  = data[indexb*7+2]-data[indexa*7+2];
  double dzi  = data[indexb*7+3]-data[indexa*7+3];
  double dvxi = data[indexb*7+4]-data[indexa*7+4];
  double dvyi = data[indexb*7+5]-data[indexa*7+5];
  double dvzi = data[indexb*7+6]-data[indexa*7+6];

  double lxi = dyi*dvzi-dzi*dvyi;
  double lyi = dzi*dvxi-dxi*dvzi;
  double lzi = dxi*dvyi-dyi*dvxi;

  double l2i = lxi*lxi + lyi*lyi + lzi*lzi;
  double li = sqrt(l2i);

  index.clear();
  index.push_back(indexa);
  index.push_back(indexb);
  vector<double> cm = get_cm_multiple(data, index);

  double dxo  = data[indexc*7+1]-cm[0*7+1];
  double dyo  = data[indexc*7+2]-cm[0*7+2];
  double dzo  = data[indexc*7+3]-cm[0*7+3];
  double dvxo = data[indexc*7+4]-cm[0*7+4];
  double dvyo = data[indexc*7+5]-cm[0*7+5];
  double dvzo = data[indexc*7+6]-cm[0*7+6];

  double lxo = dyo*dvzo-dzo*dvyo;
  double lyo = dzo*dvxo-dxo*dvzo;
  double lzo = dxo*dvyo-dyo*dvxo;

  double l2o = lxo*lxo + lyo*lyo + lzo*lzo;
  double lo = sqrt(l2o);  

  double lidotlo = lxi*lxo + lyi*lyo + lzi*lzo;

  double theta_int = acos(lidotlo/(li*lo));

  vector<double> orbelem(15);
  orbelem[0] = a1; 
  orbelem[1] = e1;
  orbelem[2] = Eb1;
  orbelem[3] = a2;
  orbelem[4] = e2;
  orbelem[5] = Eb2;
  orbelem[6] = f_stab;
  orbelem[7] = theta_int; 
  orbelem[8] = e;
  orbelem[9] = l2;
  orbelem[10] = mu;
  orbelem[11] = rcm;
  orbelem[12] = vcm;
  orbelem[13] = rcmdotvcm;
  orbelem[14] = E_ext;

  return orbelem;
}

vector<double> Products::get_binary_internal_properties(vector<double> &data, int index1, int index2) {
// a, e, EB, E_int, L2_int, mu 
  double m1  = data[index1*7];
  double x1  = data[index1*7+1];
  double y1  = data[index1*7+2];
  double z1  = data[index1*7+3];
  double vx1 = data[index1*7+4];
  double vy1 = data[index1*7+5];
  double vz1 = data[index1*7+6];

  double m2  = data[index2*7];
  double x2  = data[index2*7+1];
  double y2  = data[index2*7+2];
  double z2  = data[index2*7+3];
  double vx2 = data[index2*7+4];
  double vy2 = data[index2*7+5];
  double vz2 = data[index2*7+6];

  double mu = m1+m2;

  double x = x2-x1;
  double y = y2-y1;
  double z = z2-z1;
  double vx = vx2-vx1;
  double vy = vy2-vy1;
  double vz = vz2-vz1;  
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double v2 = vx*vx + vy*vy + vz*vz;

  double e = 0.5*v2 - mu/r;
  double lx = y*vz-z*vy;
  double ly = z*vx-x*vz;
  double lz = x*vy-y*vx;  
  double l2 = lx*lx+ly*ly+lz*lz;

  double mya = -1;
  double mye = -1;
  double myEb = 1e-100;

  double D = mu*mu + 2.0*e*l2;
  if(D >= 0 && e < 0) {
    double rp = (-mu+sqrt(D)) / (2.0*e);  
    double ra = (-mu-sqrt(D)) / (2.0*e);    

    mya = 0.5*(rp+ra);
    mye = 1.0-2.0/(ra/rp+1.0);
    myEb = 0.5*m1*m2/mya;
  }

  vector<double> orbelem(6);
  orbelem[0] = mya; 
  orbelem[1] = mye;
  orbelem[2] = myEb;
  orbelem[3] = e;
  orbelem[4] = l2;
  orbelem[5] = mu;

  return orbelem;
}

/////////////////////////////////////////////////////////////////////

void Products::detect_escapers(vector<double> &data, double Resc) {
  Nesc = 0;
  index_esc.clear();
  data_esc.clear(); 

  double Resc2 = Resc*Resc;

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    // r > Resc
    double x = data[i*7+1];
    double y = data[i*7+2];
    double z = data[i*7+3];
    double r2 = x*x + y*y + z*z;
    if(r2 > Resc2) {
      // rdotv > 0
      double vx = data[i*7+4];
      double vy = data[i*7+5];
      double vz = data[i*7+6];        
      double rdotv = x*vx + y*vy + z*vz;
      if(rdotv > 0) {
        // e >= 0;
        double e = get_energy_particle(data, i);
	if(e >= 0) {
          Nesc++;
          index_esc.push_back(i);

          vector<double> property = get_escaper_properties(data, i);
          data_esc.push_back(property);
	}
      }
    }
  }  
}

void Products::detect_binaries(vector<double> &data) {
  detect_binaries(data, this->Resc, this->fiso_bin);
}
void Products::detect_binaries(vector<double> &data, double Resc, double fiso_bin) {
  Nbin = 0;
  index_bin.clear();
  data_bin.clear();
  Nbin_esc = 0;
  index_bin_esc.clear();
  data_bin_esc.clear();

  int N = data.size()/7;
  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {

      // energy < 0
      vector<double> orbelem = get_binary_properties(data, i, j);
      double e = orbelem[3];
      if(e < 0) {

	// Isolated orbit
        double a = orbelem[0];
        vector<int> index;
        index.push_back(i);
        index.push_back(j);
        double Rnn = get_distance_nn_multiple(data, index); 
        if( Rnn/a > fiso_bin ) {
          Nbin++;
  	  index_bin.push_back(index);
          data_bin.push_back(orbelem);

          double R = orbelem[6];
          if(R > Resc) {
            double RdotV = orbelem[8];
            if(RdotV > 0) {
              double E_ext = orbelem[9];
              if(E_ext >= 0) {
                Nbin_esc++;
  	        index_bin_esc.push_back(index);
                data_bin_esc.push_back(orbelem);
              }
            }
          }         
        }
      }
    }
  } 
}

void Products::detect_2body_encounters(vector<double> &data) {
  Nbin = 0;
  index_bin.clear();
  data_bin.clear();
  Nbin_esc = 0;
  index_bin_esc.clear();
  data_bin_esc.clear();

  Nflyby = 0;
  index_flyby.clear(); 
  data_flyby.clear();  

  int N = data.size()/7;
  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {

      // energy < 0
      vector<double> orbelem = get_binary_properties(data, i, j);
      double e = orbelem[3];
      if(e < 0) {

	// Isolated orbit
        double a = orbelem[0];
        vector<int> index;
        index.push_back(i);
        index.push_back(j);
        double Rnn = get_distance_nn_multiple(data, index); 
        if( Rnn/a > fiso_bin ) {
          Nbin++;
  	  index_bin.push_back(index);
          data_bin.push_back(orbelem);

          double R = orbelem[6];
          if(R > Resc) {
            double RdotV = orbelem[8];
            if(RdotV > 0) {
              double E_ext = orbelem[9];
              if(E_ext >= 0) {
                Nbin_esc++;
  	        index_bin_esc.push_back(index);
                data_bin_esc.push_back(orbelem);
              }
            }
          }         
        }

      }
      // energy >= 0
      else {

        // approaching
        double rdotv = orbelem[12];
        if(rdotv < 0) {
          // isolated flyby
          double r = orbelem[13];
          vector<int> index;
          index.push_back(i);
          index.push_back(j);
          double Rnn = get_distance_nn_multiple(data, index); 
          if( Rnn/r > fiso_flyby ) {
            Nflyby++;
  	    index_flyby.push_back(index);
            data_flyby.push_back(orbelem);
          }          
        }

      }
    }
  } 
}

void Products::detect_3body_encounters(vector<double> &data) {
  ;
}

void Products::detect_4body_encounters(vector<double> &data) {
  ;
}

/////////////////////////////////////////////////////////////////////

void Products::process_4body(vector<double> &data) {
  detect_escapers(data, Resc);
  detect_binaries(data, Resc, fiso_bin);

  Ntri_esc = 0;
  index_tri_esc.clear();
  data_tri_esc.clear();

  if(Nesc == 0 && Nbin_esc == 0) {	// 0 = Quadruple
    outcome = 0;
    isDissolved = false;
  }
  else if(Nesc == 1 && Nbin_esc == 0) {	// 1 = Triple+Escaper
    outcome = 1;
    isDissolved = false;
 
    int index1 = 0;
    int index2 = 1;
    int index3 = 2;
    if(index_esc[0] == 0) {
      index1 = 1;
      index2 = 2;
      index3 = 3;
    }
    else if(index_esc[0] == 1) {
      index1 = 0;
      index2 = 2;
      index3 = 3;
    }
    else if(index_esc[0] == 2) {
      index1 = 0;
      index2 = 1;
      index3 = 3;   
    }

    vector<double> triple_data = get_triple_properties(data, index1, index2, index3);

    double E_int = triple_data[8];
    if(E_int < 0) {

      double f_stab = triple_data[6];    
      if(f_stab > fiso_bin) {

        vector<int> index(3);
        index[0] = index1;
        index[1] = index2;
        index[2] = index3;
        if(isIsolated_rel_multiple(data, index, fiso_bin)) {

          double Rcm = triple_data[11];
          if(Rcm > Resc) {
            double RdotV = triple_data[13];
            if(RdotV > 0) {
              double E_ext = triple_data[14];
              if(E_ext >= 0) {
                isDissolved = true;

                Ntri_esc = 1;
                index_tri_esc.push_back(index);
                data_tri_esc.push_back(triple_data);
              }
            }
          }
        }
      }
    }
  }
  else if(Nesc == 0 && Nbin_esc == 2) {	// 2 = 2 Binaries
    outcome = 2;
    isDissolved = true;
  }
  else if(Nesc == 2 && Nbin_esc == 1) {	// 3 = 2 Escapers+Binary
    outcome = 3;
    isDissolved = true;
  }
  else {				// -1 = Else
    outcome = -1;
    isDissolved = false;
  }

}

/*
void Products::process(vector<double> &data) {
  detect_escapers(data);
  detect_binaries(data);

  int N = data.size()/7;
  if(N == 3) {
    if(Nesc == 0 && Nbin_esc == 0) {		// 0 = Triple
      outcome = 0;
      isDissolved = false;
    }
    else if(Nesc == 1 && Nbin_esc == 1) {	// 1 = Binary+Escaper
      outcome = 1;
      isDissolved = true;
    }
    else {					// -1 = Else
      outcome = -1;
      isDissolved = false;
    }
  }
  else if(N == 4) {
    if(Nesc == 0 && Nbin_esc == 0) {		// 0 = Quadruple
      outcome = 0;
      isDissolved = false;
    }
    else if(Nesc == 1 && Nbin_esc == 0) {	// 1 = Triple+Escaper
      outcome = 1;
      isDissolved = false;
    }
    else if(Nesc == 0 && Nbin_esc == 2) {	// 2 = 2 Binaries
      outcome = 2;
      isDissolved = true;
    }
    else if(Nesc == 2 && Nbin_esc == 1) {	// 3 = 2 Escapers+Binary
      outcome = 3;
      isDissolved = true;
    }
    else {					// -1 = Else
      outcome = -1;
      isDissolved = false;
    }
  }  
  else if(N == 5) {
    if(Nesc == 0 && Nbin_esc == 0) {		// 0 = Quintuple
      outcome = 0;
      isDissolved = false;
    }
    else if(Nesc == 1 && Nbin_esc == 0) {	// 1 = Quadruple+Escaper
      outcome = 1;
      isDissolved = false;
    }
    else if(Nesc == 1 && Nbin_esc == 2) {	// 2 = Escaper + 2 Binaries
      outcome = 2;
      isDissolved = true;
    }
    else if(Nesc == 0 && Nbin_esc == 1) {	// 3 = Binary+Triple
      outcome = 3;
      isDissolved = false;
    }
    else if(Nesc == 2 && Nbin_esc == 0) {	// 4 = 2 Escapers+Triple
      outcome = 4;
      isDissolved = false;
    }
    else if(Nesc == 3 && Nbin_esc == 1) {	// 5 = 3 Escapers+Binary
      outcome = 5;
      isDissolved = true;
    }
    else {					// -1 = Else
      outcome = -1;
      isDissolved = false;
    }
  }
  else {
    cerr << "Dissolution not determined for this N!" << endl;
    isDissolved = false;
  }
}
*/

int Products::get_Nesc() {
  return Nesc;
}
vector<int> Products::get_index_esc() {
  return index_esc;
}
vector< vector<double> > Products::get_data_esc() {
  return data_esc;
} 

int Products::get_Nbin() {
  return Nbin;
}
vector< vector<int> > Products::get_index_bin() {
  return index_bin;
}
vector< vector<double> > Products::get_data_bin() {
  return data_bin;
} 

int Products::get_Nbin_esc() {
  return Nbin_esc;
}
vector< vector<int> > Products::get_index_bin_esc() {
  return index_bin_esc;
}
vector< vector<double> > Products::get_data_bin_esc() {
  return data_bin_esc;
} 

int Products::get_Nflyby() {
  return Nflyby;
}
vector< vector<int> > Products::get_index_flyby() {
  return index_flyby;
}
vector< vector<double> > Products::get_data_flyby() {
  return data_flyby;
} 

int Products::get_Ntri_esc() {
  return Ntri_esc;
}
vector< vector<int> > Products::get_index_tri_esc() {
  return index_tri_esc;
}
vector< vector<double> > Products::get_data_tri_esc() {
  return data_tri_esc;
} 

bool Products::get_isDissolved() {
  return isDissolved;
}
int Products::get_outcome() {
  return outcome;
}

void Products::remove_particles(vector<double> &data, vector<int> &index) {
  int N = data.size()/7;
  int M = index.size();

  vector<double> data_new;
  for(int i=0; i<N; i++) {
    bool toAdd = true;
    for(int j=0; j<M; j++) {
      if(i == index[j]) {
        toAdd = false;
        break;
      }
    }
    if(toAdd) {
      for(int j=0; j<7; j++) {
        data_new.push_back(data[i*7+j]);
      }
    }
  }

  data = data_new;
}
void Products::add_particles(vector<double> &data, vector<double> &d) {
  int M = d.size()/7;
  for(int i=0; i<M; i++) {
    for(int j=0; j<7; j++) {
      data.push_back(d[i*7+j]);
    }
  }
}
void Products::add_particles(vector<double> &data, vector< vector<double> > &d) {
  int M = d.size();
  for(int i=0; i<M; i++) {
    for(int j=0; j<d[i].size(); j++) {
      data.push_back(d[i][j]);
    }
  }
}
void Products::replace_cm_particles(vector<double> &data, vector< vector<double> > &coordinates_bin) {
  int N = data.size()/7;
  int M = coordinates_bin.size();

  for(int i=0; i<M; i++) {
    int index = N-M+i;
    for(int k=0; k<2; k++) {
      for(int j=1; j<7; j++) {
        coordinates_bin[i][k*7+j] += data[index*7+j];
      }
    }
  }

  for(int i=0; i<M; i++) {
    for(int j=0; j<7; j++) {
      data.pop_back();
    }
  }

  add_particles(data, coordinates_bin);
}







