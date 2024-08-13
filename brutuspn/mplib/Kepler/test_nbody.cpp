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

  mpreal mu  = "1";
  mpreal a   = "1";
  mpreal e   = "0";
  mpreal i   = "0" * pi/"180";
  mpreal RA  = "0" * pi/"180";
  mpreal Arp = "0" * pi/"180";
  mpreal MA  = "0" * pi/"180";

  Kepler kepler;
  mpreal EA = kepler.convert_mean_to_eccentric_anomaly(a, e, i, RA, Arp, MA);

  vector<mpreal> data = kepler.get_cartesian_coordinates(mu, a, e, i, RA, Arp, EA);

  for(int i=0; i<data.size(); i++) {
    cout << data[i] << " ";
  }
  cout << endl;

  vector<mpreal> elem = kepler.get_orbelem(mu, data[0], data[1], data[2], data[3], data[4], data[5]);

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
