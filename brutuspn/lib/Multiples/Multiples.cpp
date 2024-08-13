#include "Multiples.h"

Multiples::Multiples() {
  R_encounter = 1e100;	// radius from core to search in, to speed process up
  f_isolated = 10.0;	// isolation factor for nn2 and nn1
}

// 2-body encounters

void Multiples::make_nn_list(vector<double> &data) {
  int N = data.size()/7;
  double R2_encounter = R_encounter*R_encounter;
  nn_index.assign(N, 0);
  nn_distance.assign(N, 1e100);
  nn2_distance.assign(N, 1e100);
  for(int i=0; i<N; i++) {
    double xi = data[i*7+1];
    double yi = data[i*7+2];
    double zi = data[i*7+3];
    double r2i = xi*xi + yi*yi + zi*zi;
    if(r2i < R2_encounter) 
      for(int j=0; j<N; j++) {
        if(i != j) {
          double xj = data[j*7+1];
          double yj = data[j*7+2];
          double zj = data[j*7+3];
          double r2j = xj*xj + yj*yj + zj*zj;
          if(r2j < R2_encounter) { 
            double dx = xj-xi;
            double dy = yj-yi;
            double dz = zj-zi;
            double dr2 = dx*dx + dy*dy + dz*dz;
            if(dr2 < nn_distance[i]) {
              nn2_index[i] = nn_index[i];
              nn2_distance[i] = nn_distance[i];
              nn_index[i] = j;
              nn_distance[i] = dr2;
            } 
          }
        }
      }
    }
  }
  for(int i=0; i<N; i++) {
    nn_distance[i] = sqrt(nn_distance[i]);
    nn2_distance[i] = sqrt(nn2_distance[i]);
  }
}

vector<double> Multiples::get_binary_properties(vector<double> &data, int index1, int index2) {
// a, e, rp, ra, EB, E, L2, Rcm, Vcm, EKcm, EPcm, Ecm, RcmdotVcm 

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

  double x = x2-x1;
  double y = y2-y1;
  double z = z2-z1;
  double vx = vx2-vx1;
  double vy = vy2-vy1;
  double vz = vz2-vz1;  
  double r2 = x*x + y*y + z*z;
  double v2 = vx*vx + vy*vy + vz*vz;
  double e = 0.5*v2 - mu/sqrt(r2);
  double lx = y*vz-z*vy;
  double ly = z*vx-x*vz;
  double lz = x*vy-y*vx;  
  double l2 = lx*lx+ly*ly+lz*lz;

  double D = mu*mu + 2.0*e*l2;
  double rp = (-mu+sqrt(D)) / (2.0*e);  
  double ra = (-mu-sqrt(D)) / (2.0*e);    

  double mya = -mu/(2.0*e);
  double mye = 1.0-2.0/(ra/rp+1.0);
  double myEB = 0.5*m1*m2/mya;

  double ekcm = 0.5*mu*v2cm;
  double epcm = 0;
  int N = data.size()/7;  
  for(int i=0; i<N; i++) {
    if(i != index1 && i != index2) {
      double dx = data[i*7+1] - xcm;
      double dy = data[i*7+2] - ycm;
      double dz = data[i*7+3] - zcm;	
      double dr2 = dx*dx + dy*dy + dz*dz;
      epcm -= data[i*7]*mymu/sqrt(dr2);
    }	    
  }
  double Ecm = ekcm + epcm;

  vector<double> orbelem(13);
  orbelem[0] = mya; 
  orbelem[1] = mye;
  orbelem[2] = rp;
  orbelem[3] = ra;
  orbelem[4] = myEB;
  orbelem[5] = e;
  orbelem[6] = l2;
  orbelem[7] = rcm;
  orbelem[8] = vcm;
  orbelem[9] = ekcm;
  orbelem[10] = epcm;
  orbelem[11] = myEcm;
  orbelem[12] = rcmdotvcm;

  return orbelem;
}

void Multiples::detect_two_body_encounters(vector<double> &data) {
  int N = nn_index.size();
  for(int i=0; i<N; i++) {
    int j1 = nn_index[i];
    int j2 = nn2_index[i];
    double dr1 = nn_distance[i];
    double dr2 = nn2_distance[i];
    if(dr2/dr1 > f_isolated) {		// particle i isolated with nn
      int k1 = nn_index[j1];
      int k2 = nn2_index[j1];
      double du1 = nn_distance[j1];
      double du2 = nn_distance[j1];
      if(du2/du1 > f_isolated) {	// particle j isolated with nn
        if(i == k1) {			// eachother's nn
          vector<double> orbelem = get_binary_properties(data, i, j1);
          double e = orbelem[5];
          if(e < 0) {			// Bound encounter
            // isolated orbit?;
            // If so, isolated binary
          }
	  else {			// Unbound encounter
            // Are they approaching? rdotv < 0
            // para/hyperbolic 2-body encounter;
	  }
        }
      }
    } 
  }
}

// 3-body encounters


// 4-body encounters


// total
void Multiples::process(vector<double> &data) {
  make_nn_list(data);
  detect_two_body_encounters(data);
}






