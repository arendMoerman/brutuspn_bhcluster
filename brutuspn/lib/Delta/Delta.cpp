#include "Delta.h"

double Delta::get_phase_space_distance(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dph2 = 0;
    for(int i=0; i<N; i++) {
      double dx = data2[i*7+1] - data1[i*7+1];
      double dy = data2[i*7+2] - data1[i*7+2];
      double dz = data2[i*7+3] - data1[i*7+3];
      double dvx = data2[i*7+4] - data1[i*7+4];
      double dvy = data2[i*7+5] - data1[i*7+5];
      double dvz = data2[i*7+6] - data1[i*7+6];
      dph2 += dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
    }
    if(dph2 == 0) dph2 = 1e-100;
    double delta = 0.5*log10(dph2);
    return delta;
  }
}
double Delta::get_phase_space_distance_normalized(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dph2 = 0;
    for(int i=0; i<N; i++) {
      double dx = data2[i*7+1] - data1[i*7+1];
      double dy = data2[i*7+2] - data1[i*7+2];
      double dz = data2[i*7+3] - data1[i*7+3];
      double dvx = data2[i*7+4] - data1[i*7+4];
      double dvy = data2[i*7+5] - data1[i*7+5];
      double dvz = data2[i*7+6] - data1[i*7+6];
      dph2 += dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
    }
    dph2 /= 6.0*N;
    if(dph2 == 0) dph2 = 1e-100;
    double delta = 0.5*log10(dph2);
    return delta;
  }
}

double Delta::get_position_space_distance(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dr2 = 0;
    for(int i=0; i<N; i++) {
      double dx = data2[i*7+1] - data1[i*7+1];
      double dy = data2[i*7+2] - data1[i*7+2];
      double dz = data2[i*7+3] - data1[i*7+3];
      dr2 += dx*dx + dy*dy + dz*dz;
    }
    if(dr2 == 0) dr2 = 1e-100;
    double delta = 0.5*log10(dr2);
    return delta;
  }
}
double Delta::get_position_space_distance_normalized(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dr2 = 0;
    for(int i=0; i<N; i++) {
      double dx = data2[i*7+1] - data1[i*7+1];
      double dy = data2[i*7+2] - data1[i*7+2];
      double dz = data2[i*7+3] - data1[i*7+3];
      dr2 += dx*dx + dy*dy + dz*dz;
    }
    dr2 /= 3.0*N;
    if(dr2 == 0) dr2 = 1e-100;
    double delta = 0.5*log10(dr2);
    return delta;
  }
}

double Delta::get_velocity_space_distance(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dv2 = 0;
    for(int i=0; i<N; i++) {
      double dvx = data2[i*7+4] - data1[i*7+4];
      double dvy = data2[i*7+5] - data1[i*7+5];
      double dvz = data2[i*7+6] - data1[i*7+6];
      dv2 += dvx*dvx + dvy*dvy + dvz*dvz;
    }
    if(dv2 == 0) dv2 = 1e-100;
    double delta = 0.5*log10(dv2);
    return delta;
  }
}
double Delta::get_velocity_space_distance_normalized(vector<double> &data1, vector<double> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  if(N1 != N2) {
    cerr << "Delta: N1 != N2" << endl;
    return 1e100;
  }
  else {
    int N = N1;
    double dv2 = 0;
    for(int i=0; i<N; i++) {
      double dvx = data2[i*7+4] - data1[i*7+4];
      double dvy = data2[i*7+5] - data1[i*7+5];
      double dvz = data2[i*7+6] - data1[i*7+6];
      dv2 += dvx*dvx + dvy*dvy + dvz*dvz;
    }
    dv2 /= 3.0*N;
    if(dv2 == 0) dv2 = 1e-100;
    double delta = 0.5*log10(dv2);
    return delta;
  }
}





