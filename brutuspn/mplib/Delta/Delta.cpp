#include "Delta.h"

Delta::Delta() {
  tolerance = log10("1e-6");
}
Delta::Delta(mpreal tolerance) {
  if(tolerance < "0") this->tolerance = tolerance;
  else this->tolerance = log10(tolerance);
}

void Delta::set_tolerance(mpreal tolerance) {
  if(tolerance < "0") this->tolerance = tolerance;
  else this->tolerance = log10(tolerance);
}

mpreal Delta::get_tolerance() {
  return tolerance;
}

mpreal Delta::get_phase_space_distance(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];
    dph2 += dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
  }
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}
mpreal Delta::get_phase_space_distance_normalized(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];
    dph2 += dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
  }
  dph2 /= "6"*(mpreal)N;
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}

mpreal Delta::get_position_space_distance(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    dph2 += dx*dx + dy*dy + dz*dz;
  }
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}
mpreal Delta::get_position_space_distance_normalized(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    dph2 += dx*dx + dy*dy + dz*dz;
  }
  dph2 /= "3"*(mpreal)N;
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}

mpreal Delta::get_velocity_space_distance(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];
    dph2 += dvx*dvx + dvy*dvy + dvz*dvz;
  }
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}
mpreal Delta::get_velocity_space_distance_normalized(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  mpreal dph2 = "0";
  for(int i=0; i<N; i++) {
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];
    dph2 += dvx*dvx + dvy*dvy + dvz*dvz;
  }
  dph2 /= "3"*(mpreal)N;
  if(dph2 == "0") return dph2;

  mpreal delta = "0.5"*log10(dph2);
  return delta;
}

vector<mpreal> Delta::get_all(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  vector<mpreal> dq(6, "0");

  for(int i=0; i<N; i++) {
    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];

    dq[0] += dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
    dq[1] += dx*dx + dy*dy + dz*dz;
    dq[2] += dvx*dvx + dvy*dvy + dvz*dvz;
  }

  dq[3] = dq[0] / ("6"*(mpreal)N);
  dq[4] = dq[1] / ("3"*(mpreal)N);
  dq[5] = dq[2] / ("3"*(mpreal)N);

  for(int i=0; i<6; i++) {
    if(dq[i] != "0") dq[i] = "0.5"*log10(dq[i]);
  }

  return dq;
}

vector< vector<mpreal> > Delta::get_all_individual(vector<mpreal> &data1, vector<mpreal> &data2) {
  int N1 = data1.size()/7;
  int N2 = data2.size()/7;
  int N = N1;
  if(N2 < N1) N = N2;

  vector< vector<mpreal> > dq;

  for(int i=0; i<N; i++) {
    vector<mpreal> mydq(6, "0");

    mpreal dx = data2[i*7+1] - data1[i*7+1];
    mpreal dy = data2[i*7+2] - data1[i*7+2];
    mpreal dz = data2[i*7+3] - data1[i*7+3];
    mpreal dvx = data2[i*7+4] - data1[i*7+4];
    mpreal dvy = data2[i*7+5] - data1[i*7+5];
    mpreal dvz = data2[i*7+6] - data1[i*7+6];

    mydq[0] = dx*dx + dy*dy + dz*dz + dvx*dvx + dvy*dvy + dvz*dvz;
    mydq[1] = dx*dx + dy*dy + dz*dz;
    mydq[2] = dvx*dvx + dvy*dvy + dvz*dvz;

    mydq[3] = mydq[0] / ("6"*(mpreal)N);
    mydq[4] = mydq[1] / ("3"*(mpreal)N);
    mydq[5] = mydq[2] / ("3"*(mpreal)N);

    for(int j=0; j<6; j++) {
      if(mydq[j] != "0") mydq[j] = "0.5"*log10(mydq[j]);
    }

    dq.push_back(mydq);
  }

  return dq;
}

bool Delta::is_converged(vector<mpreal> &data1, vector<mpreal> &data2) {
  mpreal dph = get_phase_space_distance_normalized(data1, data2);
  if(dph < tolerance) return true;
  else return false;
}
mpreal Delta::get_estimate_new_tolerance(mpreal t_begin, mpreal t_end, mpreal t_div, mpreal tolerance_begin, mpreal tolerance_end) {
  if(tolerance_begin > "0") tolerance_begin = log10(tolerance_begin);
  if(tolerance_end > "0") tolerance_end = log10(tolerance_end);	
  mpreal dt = t_div - t_begin;
  mpreal dtol = tolerance_end - tolerance_begin;
  mpreal slope = dtol/dt;
  mpreal estimated_tol = (tolerance-"3") - slope*(t_end-t_begin);
  return estimated_tol;	
}





