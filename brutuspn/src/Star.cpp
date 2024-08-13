#include "Star.h"

Star::Star() {
  m = "0";
  r.fill("0");
  v.fill("0");
  a.fill("0");
  w.fill("0");
  
  vp.fill("0");
}

Star::Star(mpreal m, array<mpreal, 3> r, array<mpreal, 3> v) {
  this->m = m;
  this->r = r;
  this->v = v;
  this->w = v;
}



