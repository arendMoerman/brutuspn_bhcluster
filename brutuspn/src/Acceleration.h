#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>
#include <fstream>

#include "Star.h"
#include "mpreal.h"
using namespace mpfr;

#ifndef __Acceleration_h
#define __Acceleration_h

class Acceleration : public Star {
public:
    vector<Star> s;
    mpreal eps2;
    mpreal dt;

    mpreal vi2, vni;
    mpreal vj2, vnj;
    mpreal dvij2;
    
    mpreal Rij, apreij;
    mpreal viInvj, vijn, vijvj;
    mpreal mi, mj;
    mpreal cl, c2, c5;
    mpreal pi = mpfr::const_pi();
    
    array<mpreal, 3> barv, bara; // For diagnostics only

    Acceleration() : Star() {
        set_cl();
        set_PN();
    }
  
    int PN1p, PN1c, PN2, PN2_5, PN3, PN3_5;
    void set_PN_prod(array<mpreal, 3> &vj, array<mpreal, 3> &vi, array<mpreal, 3> &dv, array<mpreal, 3> &dr);
  
    void set_PN();
    void set_cl();

    void Acceleration_PN1_pair(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv);
    void Acceleration_PN1_cross(const Star &si, const Star &sj, array<mpreal, 3> &da, array<mpreal, 3> &dr);
    void Acceleration_PN2(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv);
    void Acceleration_PN2_5(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv);
    void Acceleration_PN3(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv);
    void Acceleration_PN3_5(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv);
    
    void get_bary(); // Not really necessary i guess
    
    void calcAcceleration_dt();
    void calcAcceleration(mpreal dt);
    
    vector<array<mpreal, 3>> getAcceleration();
    vector<array<mpreal, 3>> getJerk();
};

#endif
    
    
