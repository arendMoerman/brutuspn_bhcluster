#include <iostream>
using namespace std;

#include <vector>

#include "/usr/include/mpreal.h"
using namespace mpfr;

#ifndef __Potential_h
#define __Potential_h

class Potential {
  public:

  vector<mpreal> get_plummer_acc(mpreal M, mpreal a, mpreal x, mpreal y, mpreal z);  

};

#endif


