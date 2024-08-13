#include "Kepler.h"

mpreal Kepler::get_eccentricity(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal r2 = x*x + y*y + z*z;
  mpreal r  = sqrt(r2);
  
  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;

  mpreal ex = (vy*lz - vz*ly)/m - x/r;
  mpreal ey = (vz*lx - vx*lz)/m - y/r;
  mpreal ez = (vx*ly - vy*lx)/m - z/r;

  mpreal e2 = ex*ex + ey*ey + ez*ez;
  mpreal e  = sqrt(e2);

  return e;
}
mpreal Kepler::get_semimajor_axis(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal e = get_eccentricity(m, x, y, z, vx, vy, vz);

  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;
  
  mpreal l2 = lx*lx + ly*ly + lz*lz;

  mpreal a = l2 / (m*("1"-e*e));

  return a;   
}
mpreal Kepler::get_inclination(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;

  mpreal l2 = lx*lx + ly*ly + lz*lz;
  mpreal l = sqrt(l2);

  mpreal kx = "0";
  mpreal ky = "0";
  mpreal kz = "1";

  mpreal i = kx*lx + ky*ly + kz*lz;
  i /= l;
  i = acos(i);

  return i;
}
mpreal Kepler::get_right_ascension_of_ascending_node(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal kx = "0";
  mpreal ky = "0";
  mpreal kz = "1";

  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;

  mpreal nx = ky*lz - kz*ly;
  mpreal ny = kz*lx - kx*lz;
  mpreal nz = kx*ly - ky*lx;

  mpreal n2 = nx*nx + ny*ny + nz*nz;
  mpreal n  = sqrt(n2);

  mpreal ix = "1";
  mpreal iy = "0";
  mpreal iz = "0";

  mpreal RA = ix*nx + iy*ny + iz*nz;
  RA /= n;
  RA = acos(RA);

  mpreal jx = "0";
  mpreal jy = "1";
  mpreal jz = "0";

  mpreal pi = acos("-1");
  mpreal dir = nx*jx + ny*jy + nz*jz;
  if(dir < 0) RA = "2"*pi - RA;    

  return RA;
}
mpreal Kepler::get_argument_of_perigee(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal kx = "0";
  mpreal ky = "0";
  mpreal kz = "1";

  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;

  mpreal nx = ky*lz - kz*ly;
  mpreal ny = kz*lx - kx*lz;
  mpreal nz = kx*ly - ky*lx;

  mpreal n2 = nx*nx + ny*ny + nz*nz;
  mpreal n  = sqrt(n2);

  mpreal r2 = x*x + y*y + z*z;
  mpreal r  = sqrt(r2);

  mpreal ex = (vy*lz - vz*ly)/m - x/r;
  mpreal ey = (vz*lx - vx*lz)/m - y/r;
  mpreal ez = (vx*ly - vy*lx)/m - z/r;
  mpreal e2 = ex*ex + ey*ey + ez*ez;
  mpreal e = sqrt(e2);

  mpreal w = nx*ex + ny*ey + nz*ez;
  w /= (n*e);
  w = acos(w);

  mpreal pi = acos("-1");
  mpreal dir = ex*kx + ey*ky + ez*kz;
  if(dir < "0") w = "2"*pi - w;

  return w;
}
mpreal Kepler::get_true_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal r2 = x*x + y*y + z*z;
  mpreal r  = sqrt(r2);
  
  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;

  mpreal ex = (vy*lz - vz*ly)/m - x/r;
  mpreal ey = (vz*lx - vx*lz)/m - y/r;
  mpreal ez = (vx*ly - vy*lx)/m - z/r;

  mpreal e2 = ex*ex + ey*ey + ez*ez;
  mpreal e  = sqrt(e2);

  mpreal theta = ex*x + ey*y + ez*z;
  theta /= (e*r);
  theta = acos(theta);

  mpreal pi = acos("-1");
  mpreal rv = x*vx + y*vy + z*vz;
  if(rv < "0") theta = "2"*pi - theta;

  return theta;
}
mpreal Kepler::get_eccentric_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal theta = get_true_anomaly(m, x, y, z, vx, vy, vz);
  mpreal e = get_eccentricity(m, x, y, z, vx, vy, vz);

  mpreal E = (e + cos(theta)) / ("1" + e*cos(theta));
  E = acos(E);
  
  mpreal pi = acos("-1");
  if(theta > pi) E = "2"*pi - E;

  return E;
}
mpreal Kepler::get_mean_anomaly(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal E = get_eccentric_anomaly(m, x, y, z, vx, vy, vz);
  mpreal e = get_eccentricity(m, x, y, z, vx, vy, vz); 

  mpreal M = E - e*sin(E);

  return M;
}

mpreal Kepler::get_eccentricity(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal e = get_eccentricity(mu, x, y, z, vx, vy, vz);

  return e;
}
mpreal Kepler::get_semimajor_axis(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal a = get_semimajor_axis(mu, x, y, z, vx, vy, vz);

  return a;
}
mpreal Kepler::get_inclination(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal i = get_inclination(mu, x, y, z, vx, vy, vz);

  return i;
}
mpreal Kepler::get_right_ascension_of_ascending_node(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal RA = get_right_ascension_of_ascending_node(mu, x, y, z, vx, vy, vz);

  return RA;
}
mpreal Kepler::get_argument_of_perigee(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal w = get_argument_of_perigee(mu, x, y, z, vx, vy, vz);

  return w;
}
mpreal Kepler::get_true_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal TA = get_true_anomaly(mu, x, y, z, vx, vy, vz);

  return TA;
}
mpreal Kepler::get_eccentric_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal EA = get_eccentric_anomaly(mu, x, y, z, vx, vy, vz);

  return EA;
}
mpreal Kepler::get_mean_anomaly(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  mpreal MA = get_mean_anomaly(mu, x, y, z, vx, vy, vz);

  return MA;
}

mpreal Kepler::get_eccentricity(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal e = get_eccentricity(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return e;
}
mpreal Kepler::get_semimajor_axis(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal a = get_semimajor_axis(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return a;
}
mpreal Kepler::get_inclination(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal i = get_inclination(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return i;
}
mpreal Kepler::get_right_ascension_of_ascending_node(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal RA = get_right_ascension_of_ascending_node(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return RA;
}
mpreal Kepler::get_argument_of_perigee(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal RP = get_argument_of_perigee(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return RP;
}
mpreal Kepler::get_true_anomaly(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal TA = get_true_anomaly(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return TA;
}
mpreal Kepler::get_eccentric_anomaly(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal EA = get_eccentric_anomaly(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return EA;
}
mpreal Kepler::get_mean_anomaly(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  mpreal MA = get_mean_anomaly(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return MA;
}

vector<mpreal> Kepler::get_orbelem(mpreal &m, mpreal &x, mpreal &y, mpreal &z, mpreal &vx, mpreal &vy, mpreal &vz) {
  mpreal r2 = x*x + y*y + z*z;
  mpreal r = sqrt(r2);

  mpreal v2 = vx*vx + vy*vy + vz*vz;

  mpreal rvdot = x*vx + y*vy + z*vz;

  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;
  mpreal l2 = lx*lx + ly*ly + lz*lz;
  mpreal l = sqrt(l2);

  mpreal ix = "1";
  mpreal iy = "0";
  mpreal iz = "0";

  mpreal kx = "0";
  mpreal ky = "0";
  mpreal kz = "1";

  mpreal nx = ky*lz - kz*ly;
  mpreal ny = kz*lx - kx*lz;
  mpreal nz = kx*ly - ky*lx;
  mpreal n2 = nx*nx + ny*ny + nz*nz;
  mpreal n = sqrt(n2);

  mpreal C1 = (v2-m/r)/m;
  mpreal C2 = rvdot/m;
  mpreal ex = C1*x - C2*vx;
  mpreal ey = C1*y - C2*vy;
  mpreal ez = C1*z - C2*vz;
  mpreal e2 = ex*ex + ey*ey + ez*ez;
  mpreal e = sqrt(e2);

  mpreal E = v2/"2" - m/r;

  mpreal a = "0";
  mpreal p = "0";
  if(e < "1") {
    a = -m/("2"*E); 
    p = a*("1"-e*e);
  }
  else {
    a = "1e100";
    p = l2/m;
  }

  mpreal i = acos(lz/l);
  mpreal RA = "0";
  mpreal Arp = "0"; // How to determine angle of pericenter direction for i=0?
  mpreal TA = "0";
  if(n == "0") {
    TA  = acos( (a*(1.-e2)/r - 1.) / e );
  }
  else {
    RA = ix*nx + iy*ny + iz*nz;
    RA /= n;
    RA = acos(RA);

    mpreal ndote = nx*ex + ny*ey + nz*ez;
    Arp = acos(ndote / (n*e));

    mpreal ndotr = nx*x + ny*y + nz*z;
    TA = acos(ndotr / (n*r));  
  }

  if(rvdot < "0") {
    TA = 2.*acos("-1") - TA;
  }

  TA = get_true_anomaly(m, x, y, z, vx, vy, vz);

  mpreal MA = convert_true_to_mean_anomaly(a, e, i, RA, Arp, TA);

  //cerr << a << " " << e << " " << i << " " << RA << " " << Arp << " " << MA << endl;

  vector<mpreal> elem(6);
  elem[0] = a;
  elem[1] = e;
  elem[2] = i;
  elem[3] = RA;
  elem[4] = Arp;
  elem[5] = MA;

  return elem;

/*  vector<mpreal> elem(8);

  mpreal r2 = x*x + y*y + z*z;
  mpreal r  = sqrt(r2);
  
  mpreal lx = y*vz - z*vy;
  mpreal ly = z*vx - x*vz;
  mpreal lz = x*vy - y*vx;
  mpreal l2 = lx*lx + ly*ly + lz*lz;
  mpreal l = sqrt(l2);

  mpreal ex = (vy*lz - vz*ly)/m - x/r;
  mpreal ey = (vz*lx - vx*lz)/m - y/r;
  mpreal ez = (vx*ly - vy*lx)/m - z/r;

  mpreal e2 = ex*ex + ey*ey + ez*ez;
  mpreal e  = sqrt(e2);
  elem[1] = e;

  mpreal a = l2 / (m*("1"-e*e));
  elem[0] = a;

  mpreal kx = "0";
  mpreal ky = "0";
  mpreal kz = "1";

  mpreal i = kx*lx + ky*ly + kz*lz;
  i /= l;
  i = acos(i);
  elem[2] = i;

  mpreal nx = ky*lz - kz*ly;
  mpreal ny = kz*lx - kx*lz;
  mpreal nz = kx*ly - ky*lx;
  mpreal n2 = nx*nx + ny*ny + nz*nz;
  mpreal n  = sqrt(n2);

  mpreal ix = "1";
  mpreal iy = "0";
  mpreal iz = "0";

  mpreal RA = ix*nx + iy*ny + iz*nz;
  RA /= n;
  RA = acos(RA);

  mpreal jx = "0";
  mpreal jy = "1";
  mpreal jz = "0";

  mpreal pi = acos("-1");
  mpreal dir = nx*jx + ny*jy + nz*jz;
  if(dir < 0) RA = "2"*pi - RA;   
  elem[3] = RA;

  mpreal w = nx*ex + ny*ey + nz*ez;
  w /= (n*e);
  w = acos(w);

  dir = ex*kx + ey*ky + ez*kz;
  if(dir < "0") w = "2"*pi - w;  
  elem[4] = w;
  
  mpreal theta = ex*x + ey*y + ez*z;
  theta /= (e*r);
  theta = acos(theta);

  mpreal rv = x*vx + y*vy + z*vz;
  if(rv < "0") theta = "2"*pi - theta;  
  elem[5] = theta;

  mpreal E = (e + cos(theta)) / ("1" + e*cos(theta));
  E = acos(E);
  
  if(theta > pi) E = "2"*pi - E;  
  elem[6] = E;

  mpreal M = E - e*sin(E);  
  elem[7] = M;

  return elem;*/
}
vector<mpreal> Kepler::get_orbelem(mpreal &m1, mpreal &x1, mpreal &y1, mpreal &z1, mpreal &vx1, mpreal &vy1, mpreal &vz1, mpreal &m2, mpreal &x2, mpreal &y2, mpreal &z2, mpreal &vx2, mpreal &vy2, mpreal &vz2) {
  mpreal mu = m1 + m2;
  mpreal x = x2-x1;
  mpreal y = y2-y1;
  mpreal z = z2-z1;
  mpreal vx = vx2-vx1;
  mpreal vy = vy2-vy1;
  mpreal vz = vz2-vz1;

  vector<mpreal> elem = get_orbelem(mu, x, y, z, vx, vy, vz);

  return elem;
}
vector<mpreal> Kepler::get_orbelem(vector<mpreal> &data, int index1, int index2) {
  mpreal m1  = data[index1*7+0];
  mpreal x1  = data[index1*7+1];
  mpreal y1  = data[index1*7+2];
  mpreal z1  = data[index1*7+3];
  mpreal vx1 = data[index1*7+4];
  mpreal vy1 = data[index1*7+5];
  mpreal vz1 = data[index1*7+6];

  mpreal m2  = data[index2*7+0];
  mpreal x2  = data[index2*7+1];
  mpreal y2  = data[index2*7+2];
  mpreal z2  = data[index2*7+3];
  mpreal vx2 = data[index2*7+4];
  mpreal vy2 = data[index2*7+5];
  mpreal vz2 = data[index2*7+6];

  vector<mpreal> orbelem = get_orbelem(m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2);

  return orbelem;
}

mpreal Kepler::convert_eccentric_to_true_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal TA = "2" * atan2(sqrt("1"+e)*sin(A/"2"), sqrt("1"-e)*cos(A/"2"));
  return TA;
}
mpreal Kepler::convert_true_to_eccentric_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal EA = atan2(sqrt("1"-e*e)*sin(A), "e"+cos(A));
  return EA;
}
mpreal Kepler::convert_eccentric_to_mean_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal MA = A - e*sin(A);
  return MA;
}
mpreal Kepler::convert_mean_to_eccentric_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal delta = "1e-15";

  mpreal EA = A;

  mpreal diff = "1";
  if(e > "0.4") {
    while(diff > delta) {
      mpreal EA0 = EA;
      EA = EA0 + (A + e*sin(EA0) - EA0) / ("1" - e*cos(EA0));
      diff = abs(EA-EA0);
    }
  }
  else {
    while(diff > delta) {
      mpreal EA0 = EA;
      EA = A + e*sin(EA0);
      diff = abs(EA-EA0);
    }
  }

  return EA;  
}
mpreal Kepler::convert_mean_to_true_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal EA = convert_mean_to_eccentric_anomaly(a, e, i, RA, Arp, A);
  mpreal TA = convert_eccentric_to_true_anomaly(a, e, i, RA, Arp, EA);
  return TA;
}
mpreal Kepler::convert_true_to_mean_anomaly(mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &A) {
  mpreal EA = convert_true_to_eccentric_anomaly(a, e, i, RA, Arp, A);
  mpreal MA = convert_eccentric_to_mean_anomaly(a, e, i, RA, Arp, EA);
  return MA;
}

vector<mpreal> Kepler::get_cartesian_coordinates(mpreal &mu, mpreal &a, mpreal &e, mpreal &i, mpreal &RA, mpreal &Arp, mpreal &E) {
  mpreal px = cos(Arp) * cos(RA) - sin(Arp) * cos(i) * sin(RA);
  mpreal py = cos(Arp) * sin(RA) + sin(Arp) * cos(i) * cos(RA);
  mpreal pz = sin(i) * sin(Arp);

  mpreal qx = -sin(Arp) * cos(RA) - cos(Arp) * cos(i) * sin(RA);
  mpreal qy = -sin(Arp) * sin(RA) + cos(Arp) * cos(i) * cos(RA);
  mpreal qz = sin(i) * cos(Arp);

  mpreal C1 = a*(cos(E)-e);
  mpreal C2 = a*sqrt("1"-e*e)*sin(E);
  mpreal x = C1*px + C2*qx;
  mpreal y = C1*py + C2*qy;
  mpreal z = C1*pz + C2*qz;  

  mpreal Edot = sqrt(mu/(a*a*a)) / ("1"-e*cos(E));
  mpreal C3 = -a*sin(E)*Edot;
  mpreal C4 = a*sqrt("1"-e*e)*cos(E)*Edot;
  mpreal vx = C3*px + C4*qx;
  mpreal vy = C3*py + C4*qy;
  mpreal vz = C3*pz + C4*qz;

  vector<mpreal> data(6);
  data[0] = x;  
  data[1] = y;  
  data[2] = z;  
  data[3] = vx;  
  data[4] = vy;  
  data[5] = vz;

  return data;  
}



