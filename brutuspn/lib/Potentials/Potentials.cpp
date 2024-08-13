#include "Potentials.h"

Potentials::Potentials() {
  potential_type = "none";
  tolerance = 1e-14;
  Mgal = 1;
}

void Potentials::set_potential_type(string type, vector<double> parm) {
  potential_type = type;
  par = parm;
}
void Potentials::set_tolerance(double tol) {
  tolerance = tol;
}

double Potentials::get_point_mass_potential_at(double m, double x, double y, double z) {
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double phi = -Mgal*m/r;
  return phi;
}
double Potentials::get_homo_sphere_potential_at(double m, double x, double y, double z) {
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double pi = acos(-1.0);
  double rho = Mgal/(4./3*pi*pow(par[0], 3));
  double phi;
  if(r < par[0]) phi = -2.*pi*rho*(par[0]*par[0]-r2/3)*m;
  else phi = -4.*pi*rho*par[0]*par[0]*par[0]/(3.*r)*m;
  return phi;
}
double Potentials::get_isochrone_potential_at(double m, double x, double y, double z) {
  double r2 = x*x + y*y + z*z;
  double r = sqrt(r2);
  double phi = -Mgal*m/(par[0] + sqrt(r*r + par[0]*par[0]));
  return phi;
}

vector<double> Potentials::get_acceleration_at(double m, double x, double y, double z) {
  vector<double> acc(3, 0);

  double phi0 = 0;
  if(potential_type == "point_mass") phi0 = get_point_mass_potential_at(m, x, y, z);
  else if(potential_type == "homo_sphere") phi0 = get_homo_sphere_potential_at(m, x, y, z);
  else if(potential_type == "isochrone") phi0 = get_isochrone_potential_at(m, x, y, z);

  double dq0 = 0.125;
  double dq;
  int n = 1;

  vector<double> dqs;
  vector< vector<double> > as;
  vector<double> a;
  vector<double> a_exp0;
  vector<double> a_exp1;
  bool flag = false;

  // n=1
  dq = dq0/n;
  a = iterate(phi0, dq, m, x, y, z);
  dqs.push_back(dq);  
  as.push_back(a);

  extrapol(dqs, as, a_exp0);

  // n=2
  n = 2;
  dq = dq0/n;
  a = iterate(phi0, dq, m, x, y, z);
  dqs.push_back(dq);  
  as.push_back(a);

  extrapol(dqs, as, a_exp1);

  flag = error_control(a_exp0, a_exp1);

  if(!flag) {
    n += 2;
    dq = dq0/n;
    a = iterate(phi0, dq, m, x, y, z);
    dqs.push_back(dq);  
    as.push_back(a);

    a_exp0 = a_exp1;
    extrapol(dqs, as, a_exp1);

    flag = error_control(a_exp0, a_exp1);
  }

  acc = a_exp1;

  return acc;
}
vector<double> Potentials::iterate(double phi0, double dq, double m, double x, double y, double z) {
  vector<double> acc(3, 0);

  double phix = 0;
  double phiy = 0;
  double phiz = 0;
  if(potential_type == "point_mass") {
    phix = get_point_mass_potential_at(m, x+dq, y, z);
    phiy = get_point_mass_potential_at(m, x, y+dq, z);
    phiz = get_point_mass_potential_at(m, x, y, z+dq);
  }
  else if(potential_type == "homo_sphere") {
    phix = get_homo_sphere_potential_at(m, x+dq, y, z);
    phiy = get_homo_sphere_potential_at(m, x, y+dq, z);
    phiz = get_homo_sphere_potential_at(m, x, y, z+dq);
  }
  else if(potential_type == "isochrone") {
    phix = get_isochrone_potential_at(m, x+dq, y, z);
    phiy = get_isochrone_potential_at(m, x, y+dq, z);
    phiz = get_isochrone_potential_at(m, x, y, z+dq);
  }

  double dphi_dx = (phix-phi0)/dq;
  double dphi_dy = (phiy-phi0)/dq;
  double dphi_dz = (phiz-phi0)/dq;

  acc[0] = -dphi_dx/m;
  acc[1] = -dphi_dy/m;
  acc[2] = -dphi_dz/m;

  return acc;
}
void Potentials::extrapol(vector<double> &dqs, vector< vector<double> > &as, vector<double> &a_exp) {
  int M = dqs.size();
  vector<double> ax_sample(M), ay_sample(M), az_sample(M);
  a_exp.resize(3);
  for(int j=0; j<M; j++) {
    ax_sample[j]  = as[j][0];
    ay_sample[j]  = as[j][1];
    az_sample[j]  = as[j][2];
  }
  a_exp[0]  = extrapolate(dqs, ax_sample, 0.0);
  a_exp[1]  = extrapolate(dqs, ay_sample, 0.0);
  a_exp[2]  = extrapolate(dqs, az_sample, 0.0);
}
double Potentials::extrapolate(vector<double> x, vector<double> y, double x0) {
  int N = x.size();
  if(N == 1) {
    return y[0];
  }
  else {
    for(int i=1; i<N; i++) {
      for(int j=0; j<N-i; j++) {
	y[j] = ( (x0-x[j+i])*y[j] + (x[j]-x0)*y[j+1] ) / ( x[j]-x[j+i] );
      }
    }
    return y[0];
  }  
}
bool Potentials::error_control(vector<double> &a1, vector<double> &a2) {
  bool flag = true;
  if( fabs(a1[0]-a2[0]) > tolerance ) flag = false;
  if( fabs(a1[1]-a2[1]) > tolerance ) flag = false;
  if( fabs(a1[2]-a2[2]) > tolerance ) flag = false;
  return flag;
}






