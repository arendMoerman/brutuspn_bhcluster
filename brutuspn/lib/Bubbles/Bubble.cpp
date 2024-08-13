#include "Bubble.h"

Bubble::Bubble() {
  N = 0;			
  data.clear();		
  data_cm.clear();	
  R = 0;	
  f_radius = 1.0;		
  dt = 0;			
}
Bubble::Bubble(double &f_radius) {
  N = 0;			
  data.clear();		
  data_cm.clear();	
  R = 0;	
  this->f_radius = f_radius;		
  dt = 0;
}
Bubble::Bubble(vector<double> &raw_data) {
  process_raw_data(raw_data);
  f_radius = 1.0;
  dt = 0;
} 
Bubble::Bubble(double &f_radius, vector<double> &raw_data) {
  process_raw_data(raw_data);
  this->f_radius = f_radius;
  dt = 0;
}  
void Bubble::set_N(int &N) {
  this->N = N;
}
void Bubble::set_data(vector<double> &data) {
  this->data = data;
}
void Bubble::set_data_cm(vector<double> &data_cm) {
  this->data_cm = data_cm;
}
void Bubble::set_R(double &R) {
  this->R = R;
}
void Bubble::set_f_radius(double &f_radius) {
  this->f_radius = f_radius;
}
void Bubble::set_dt(double &dt) {
  this->dt = dt;
}
int Bubble::get_N() {
  return N;
}
vector<double> Bubble::get_data() {
  return data;
}
vector<double> Bubble::get_data_cm() {
  return data_cm;
}
double Bubble::get_R() {
  return R;
}
double Bubble::get_f_radius() {
  return f_radius;
}
double Bubble::get_dt() {
  return dt;
}
void Bubble::process_raw_data(vector<double> &raw_data) {
  N = raw_data.size()/7;

  data_cm.assign(7, 0);
  for(int i=0; i<N; i++) {
    data_cm[0] += raw_data[i*7];
    for(int j=1; j<7; j++) {
      data_cm[j] += raw_data[i*7]*raw_data[i*7+j];
    }
  }
  for(int i=1; i<7; i++) {
    data_cm[i] /= data_cm[0];
  }

  data.assign(7*N, 0);		
  for(int i=0; i<N; i++) {
    data[i*7] = raw_data[i*7];
    for(int j=1; j<7; j++) {
      data[i*7+j] = raw_data[i*7+j]-data_cm[j];
    }
  }	

  R = 0;
  for(int i=0; i<N; i++) {
    double myR2 = 0;
    for(int j=1; j<3; j++) {
      myR2 += data[i*7+j]*data[i*7+j];
    }
    if(myR2 > R) {
      R = myR2;
    }
  }		
  R = f_radius*sqrt(R);	
}
vector<double> Bubble::get_raw_data() {
  vector<double> raw_data(7*N, 0);
  for(int i=0; i<N; i++) {
    raw_data[i*7] = data[i*7];
    for(int j=1; j<7; i++) {
      raw_data[i*7+j] = data[i*7+j]+data_cm[j];
    }
  }
  return raw_data;
}



