#include <cstdlib>
#include <cmath>
#include "mtrand.h"

#ifndef __Random_h
#define __Random_h

class Random {

  int seed, pivot;

  MTRand_int32 irand;
  MTRand drand;

  public:

  Random();
  Random(int seed);
  Random(int seed, int pivot);

  void set_seed(int seed);
  void set_pivot(int pivot);

  int uniform(int min, int max);

  float uniform(float min, float max);

  double random();
  double uniform(double min, double max);
  double gaussian(double mu, double sigma);
  double gaussian();

};

#endif


