#include "Random.h"

Random::Random() {
  seed = 0;
  pivot = 0;
  irand.seed(seed);
  drand.seed(seed);
}
Random::Random(int seed) {
  this->seed = seed;
  pivot = 0;
  irand.seed(seed);
  drand.seed(seed);
}
Random::Random(int seed, int pivot) {
  this->seed = seed;
  this->pivot = pivot;
  irand.seed(seed);
  drand.seed(seed);
  for(int i=0; i<pivot; i++) {
    int dummy = uniform(1, 7);
    double fdummy = uniform(0.0, 10.0);
  }
}

void Random::set_seed(int seed) {
  irand.seed(seed);
  drand.seed(seed);
}
void Random::set_pivot(int pivot) {
  for(int i=0; i<pivot; i++) {
    int dummy = uniform(1, 7);
    double fdummy = uniform(0.0, 10.0);
  }
}

int Random::uniform(int min, int max) {
  return irand()%(max-min)+min;
}

float Random::uniform(float min, float max) {
  return (float)drand()*(max-min)+min;
}

double Random::random() {
  return drand();
}
double Random::uniform(double min, double max) {
  return drand()*(max-min)+min;
}

double Random::gaussian(double mu, double sigma) {
  double z = 0;

  double pi = acos(-1.0);
  bool accept = false;
  while(!accept) {
    double x = drand()*20 - 10;
    double u = exp(-0.5*x*x) / sqrt(2*pi);
    double p = drand()*0.4;
    if(p < u) {
      z = x;      
      accept = true;
    }
  }

  return z*sigma+mu;
}
double Random::gaussian() {
  return gaussian(0, 1);
}






