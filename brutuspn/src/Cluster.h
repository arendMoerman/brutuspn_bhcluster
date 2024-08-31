#include <iostream>
using namespace std;

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>

#include "Star.h"
#include "Acceleration.h"

#ifndef __Cluster_h
#define __Cluster_h

struct mergerOut {
    bool merged;
    int idx0;
    int idx1;
};

class Cluster : public Acceleration { 
public:
    mpreal time, dt_last;
    Cluster() : Acceleration() {}
    
    Cluster(vector<double> data, mpreal r_merge);
    Cluster(vector<mpreal> data, mpreal r_merge);

    mpreal rmerge;

    int merge_idx0, merge_idx1;

    vector<double> get_data_double();
    vector<mpreal> get_data();   
    
    array<mpreal, 3> a0;
    
    void setw0();

    void updatePositions(mpreal dt);
    void updateAuxiliary(mpreal dt);
    void updateVelocities(mpreal dt);
    void updateAuxiliary2(mpreal dt);

    void step(mpreal &dt);
  
    mpreal e1(const Star &si, int N);
    mpreal e2(int N);
    mpreal e3(int N);
    mpreal energies();
    vector<mpreal> get_Ener();
    
    void collisionDetection(mergerOut &out);
};

#endif


