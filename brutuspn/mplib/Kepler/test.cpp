#include "iostream"
using namespace std;

#include <cstdlib>

#include "mpreal.h"
using namespace mpfr;

#include "Kepler.h"

int main(int argc, char* argv[]) {

  int numBits = atoi(argv[1]);
  mpreal::set_default_prec(numBits);
  cout.precision(numBits/4);

  mpreal pi = acos("-1");

  mpreal m1  = "1";
  mpreal m2  = "0.0009546";
  mpreal M = m1+m2;

  mpreal a   = "5.203079444926453";
  mpreal e   = "4.884278853265638E-02";
  mpreal i   = (mpreal)"1.303452576455068" / "180"*pi;
  mpreal RA  = (mpreal)"1.005209569701123E+02" / "180"*pi;
  mpreal Arp = (mpreal)"2.740412149815173E+02" / "180"*pi;
  mpreal MA  = (mpreal)"5.769452833039548E+01" / "180"*pi;

  vector<mpreal> cm(6);
  cm[0] = "-1.037225476964814E-03";
  cm[1] = "-2.500308298077913E-03";
  cm[2] = "-4.736948118279435E-05";
  cm[3] = "6.212585408179171E-06";
  cm[4] = "-9.398699675021308E-07";
  cm[5] = "-1.348510476683671E-07";

  vector<mpreal> nasa(6);
  nasa[0] = "1.128258136939083E+00";
  nasa[1] = "4.946562536985701E+00";
  nasa[2] = "-4.587250036267861E-02";
  nasa[3] = "-7.449327901882835E-03";
  nasa[4] = "2.037693791638855E-03";
  nasa[5] = "1.581827242429274E-04";

  Kepler kepler;

  mpreal E = kepler.convert_mean_to_eccentric_anomaly(a, e, i, RA, Arp, MA);

  mpreal mu = "0.00029619468";

  vector<mpreal> data = kepler.get_cartesian_coordinates(mu, a, e, i, RA, Arp, E);
  for(int i=0; i<data.size(); i++) {
    data[i] += cm[i];
  }

  for(int i=0; i<data.size(); i++) {
    cout << data[i] << endl;
    cout << nasa[i] << endl;
  }
  cout << endl;

  vector<mpreal> elem = kepler.get_orbelem(M, data[0], data[1], data[2], data[3], data[4], data[5]);

  cout << a << endl;
  cout << elem[0] << endl;
  cout << e << endl;
  cout << elem[1] << endl;
  cout << i << endl;
  cout << elem[2] << endl;
  cout << RA << endl;
  cout << elem[3] << endl;
  cout << Arp << endl;
  cout << elem[4] << endl;
  cout << MA << endl;
  cout << elem[5] << endl;

  return 0;
}
