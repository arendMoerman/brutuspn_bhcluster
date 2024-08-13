#include "Potential.h"

vector<mpreal> Potential::get_plummer_acc(mpreal M, mpreal a, mpreal x, mpreal y, mpreal z) {
  mpreal r2 = x*x + y*y + z*z;
  mpreal r = sqrt(r2);

  mpreal Menc = M * r2*r / pow(r2 + a*a, "1.5");

  mpreal acc = "0";
  mpreal ax = "0";
  mpreal ay = "0";
  mpreal az = "0";

  if(r2 > "0") {
    acc = Menc / r2;
    ax = -x/r * acc;
    ay = -y/r * acc;
    az = -z/r * acc;
  }

  vector<mpreal> accs(3);
  accs[0] = ax;
  accs[1] = ay;
  accs[2] = az;

  return accs;
}  







