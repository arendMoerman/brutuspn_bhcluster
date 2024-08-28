#include "Cluster.h"

Cluster::Cluster(vector<double> data, mpreal r_merge) {
    int N = data.size()/7;
    s.resize(N);
    mpreal m;
    array<mpreal, 3> r, v;
    for(int i=0; i<N; i++) {
        m    = (mpreal)data[i*7+0];
        r[0] = (mpreal)data[i*7+1];
        r[1] = (mpreal)data[i*7+2];
        r[2] = (mpreal)data[i*7+3];
        v[0] = (mpreal)data[i*7+4];
        v[1] = (mpreal)data[i*7+5];
        v[2] = (mpreal)data[i*7+6];
        s[i] = Star(m, r, v);
    }
    this->time = 0;
    this->rmerge = r_merge;
}
Cluster::Cluster(vector<mpreal> data, mpreal r_merge) {
    int N = data.size()/7;
    s.resize(N);
    mpreal m;
    array<mpreal, 3> r, v;
    for(int i=0; i<N; i++) {
        m    = data[i*7+0];
        r[0] = data[i*7+1];
        r[1] = data[i*7+2];
        r[2] = data[i*7+3];
        v[0] = data[i*7+4];
        v[1] = data[i*7+5];
        v[2] = data[i*7+6];
        s[i] = Star(m, r, v);
    }
    this->time = 0;
    this->rmerge = r_merge;
}

// This function sets the initial values for auxiliary vector w
void Cluster::setw0() {
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->w = si->v;
    }
}

// These update functions advance a step through the algorithm
void Cluster::updatePositions(mpreal dt) {
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        for(int k=0; k<3; k++) {
            si->r[k] += si->v[k]*dt/"2";
        }
    }
}

void Cluster::updateAuxiliary(mpreal dt) {
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        for(int k=0; k<3; k++) {
            si->w[k] += si->a[k]*dt/"2";
        }
        si->vp = si->v; //save current velocities
        si->v  = si->w; //pass w as new velocities for acceleration
    }
}

void Cluster::updateAuxiliary2(mpreal dt) {
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        for(int k=0; k<3; k++) {
            si->w[k] += si->a[k]*dt/"2";
        }
    }
}

void Cluster::updateVelocities(mpreal dt) {
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        for(int k=0; k<3; k++) {
            si->v[k] = si->vp[k] + si->a[k]*dt;
            
            if (si->v[k] > cl) {
                cout << "Error: v/cl ="
                     << si->v[k]/cl
                     << " > 1! Nice try guy ;)"
                     << endl;
                exit(1);
            }
        }
    }
}

// Step the cluster using Auxiliary Vector Algorithm by Mikkola & Hellstrom
void Cluster::step(mpreal &dt) {
    updatePositions(dt);
    calcAcceleration(dt);
    updateAuxiliary(dt);
    calcAcceleration(dt);
    updateVelocities(dt);
    calcAcceleration(dt);
    updateAuxiliary2(dt);
    updatePositions(dt);
}

// Here, the PN contributions up to 3PN are calculated.
// Note that the 2.5 and 3.5 terms have no conserved energy
mpreal Cluster::e1(const Star &si, int N) {
    mpreal E1cross = "0";
    mpreal E1;
    string Ns = to_string(N);
    
    array<mpreal, 3> ik;
    mpreal mk;
    mpreal ik2;
    mpreal Rik;
    
    for (vector<Star>::iterator sk = s.begin(); sk != s.end(); ++sk) {
        if (*sk == si) continue;
        
        ik2 = "0";
        for(int k=0; k<3; k++) {
            ik[k] = si.r[k] - sk->r[k];
            ik2  += ik[k]*ik[k];
        }
        mk = sk->m;
        Rik = sqrt(ik2 + eps2);
        
        E1cross += mk/Rik;
    }
    E1 = c2*mi*("0.375"*vi2*vi2/(Ns-"1") + "1.5"*vi2*mj/Rij + "0.5"*E1cross*mj/Rij -
        "0.25"*mj/Rij*("7"*viInvj + vni*vnj));
    return E1;
}
  
mpreal Cluster::e2(int N) {
    mpreal E2;
    string Ns = to_string(N);
    
    E2 = c2*c2*("0.3125"*mi*vi2*vi2*vi2/(Ns-"1") - "0.5"*mi*mi*mi*mj/(Rij*Rij*Rij) - 
        "2.375"*mi*mi*mj*mj/(Rij*Rij*Rij) +
        "0.5"*mi*mi*mj/(Rij*Rij)*("-3"*vi2 + "3.5"*vj2 + "14.5"*vni*vni - "6.5"*vni*vnj + vnj*vnj) +
        "0.25"*mi*mj/Rij*("1.5"*vni*vni*vni*vnj + "0.75"*vni*vni*vnj*vnj - "4.5"*vni*vnj*vi2 -
            "6.5"*vnj*vnj*vi2 + "10.5"*vi2*vi2 + "6.5"*vni*vni*viInvj +
            "3"*vni*vnj*viInvj - "27.5"*vi2*viInvj + "8.5"*viInvj*viInvj + "7.75"*vi2*vj2));
    return E2;
}

mpreal Cluster::e3(int N) {
    mpreal E3;
    string Ns = to_string(N);
    
    E3 = c2*c2*c2*("0.2734375"*mi*vi2*vi2*vi2*vi2/(Ns-"1") + 
        ("0.375"*mi + "469"*mj/"18")*mi*mi*mi*mj/(Rij*Rij*Rij*Rij) +
        "0.5"*mi*mi*mj*mj/(Rij*Rij*Rij)*("547"*vni*vni/"6" - "3115"*vni*vnj/"24" - pi*pi*"3.84375"*vni*vijn -
            "575"*vi2/"9" + pi*pi*"1.28125"*vijvj + "4429"*viInvj/"72") +
        "0.5"*mi*mi*mi*mj/(Rij*Rij*Rij)*("-109.25"*vni*vni + "79.25"*vni*vnj + "3"*vnj*vnj + "301"*vi2/"12" -
            "337"*viInvj/"12" + "2.5"*vj2) +
        mi*mj/Rij*("-0.3125"*vni*vni*vni*vni*vni*vnj - "0.3125"*vni*vni*vni*vni*vnj*vnj - 
            "0.15625"*vni*vni*vni*vnj*vnj*vnj +
            "1.1875"*vni*vni*vni*vnj*vi2 + "0.9375"*vni*vni*vnj*vnj*vi2 +
            "0.75"*vni*vnj*vnj*vnj*vi2 + "1.1875"*vnj*vnj*vnj*vnj*vi2 - "1.3125"*vni*vnj*vi2*vi2 -
            "2"*vnj*vnj*vi2*vi2 + "3.4375"*vi2*vi2*vi2 - "1.1875"*vni*vni*vni*vni*viInvj -
            vni*vni*vni*vnj*viInvj -
            "0.46875"*vni*vni*vnj*vnj*viInvj + "2.8125"*vni*vni*vi2*viInvj +
            "1.25"*vni*vnj*vi2*viInvj + "2.75"*vnj*vnj*vi2*viInvj - "8.6875"*vi2*vi2*viInvj -
            "0.75"*vni*vni*viInvj*viInvj + "0.3125"*vni*vnj*viInvj*viInvj + "5.125"*vi2*viInvj*viInvj +
            "0.0625"*viInvj*viInvj*viInvj - "2.8125"*vni*vni*vi2*vj2 - "0.71875"*vni*vnj*vi2*vj2 + "4.9375"*vi2*vi2*vj2 -
            "5.03125"*vi2*vj2*viInvj) +
        mi*mi*mj/(Rij*Rij)*("-6.125"*vni*vni*vni*vni + "9.375"*vni*vni*vni*vnj - "23.375"*vni*vni*vnj*vnj + 
            "5.5"*vi2*vi2 +
            "247"*vni*vnj*vnj*vnj/"24" + "6.125"*vni*vni*vi2 + "10.125"*vni*vnj*vi2 -
            "5.25"*vnj*vnj*vi2 - "7.5"*vni*vni*viInvj - "1.5"*vni*vnj*viInvj +
            "5.25"*vnj*vnj*viInvj - "27"*vi2*viInvj + "27.5"*viInvj*viInvj + "12.25"*vni*vni*vj2 -
            "13.5"*vni*vnj*vj2 + "0.75"*vnj*vnj*vj2 + "13.75"*vi2*vj2 - "28"*viInvj*vj2 + "8.4375"*vj2*vj2));
    return E3;
}

// Function that calculates the total energy in the cluster
mpreal Cluster::energies() {
    int N = s.size();
    
    mpreal E0 = "0";
    mpreal E1 = "0";
    mpreal E2 = "0";
    mpreal E3 = "0";
  
    mpreal Ekini = "0";
    mpreal E     = "0";
    
    mpreal dr2   = "0";
    
    array<mpreal, 3> vi;
    array<mpreal, 3> dr;
    array<mpreal, 3> dv;
    array<mpreal, 3> da;

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->ener = "0";
        mpreal Epotj = "0";
        vi2 = "0";
        array<mpreal, 3> vj;
        
        for(int k=0; k<3; k++) {
            vi[k] = si->v[k];
            vi2  += vi[k]*vi[k];
        }
        mi = si->m;
        
        Ekini = "0.5"*mi*vi2;
        si->ener += Ekini;
    
        for (vector<Star>::iterator sj = s.begin(); sj != s.end(); ++sj) {
            if(sj == si) continue;
            
            dr2 = "0";
            for(int k=0; k<3; k++) {
                vj[k] = sj->v[k];
                dr[k] = si->r[k] - sj->r[k];
                dv[k] = si->v[k] - sj->v[k];
                dr2  += dr[k]*dr[k];
            }
            mj = sj->m;
            Rij = sqrt(dr2 + eps2);
                        
            set_PN_prod(vj, vi, dv, dr);
            
            E0 = -mi*mj/("2"*Rij);
            
            if (PN1p)   E1 = e1(*si, N);
            if (PN2)    E2 = e2(N);
            if (PN3)    E3 = e3(N);
            
            Epotj = (E0 + E1 + E2 + E3);
            
            E += Epotj;
            si->ener += Epotj;
        }
        E += Ekini;
    }
    return E;
}

// Function that returns a vector with the individual energies
vector<mpreal> Cluster::get_Ener() {
    vector<mpreal> Ener;
    
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        Ener.push_back(si->ener);
    }
    return Ener;
}

// Function to detect a collision for length scale = 1 parsec!!
void Cluster::collisionDetection(mergerOut &out) {
    out.merged = false;
    out.idx0 = -1;
    out.idx1 = -1;
    
    mpreal dx, dy, dz;
    mpreal dr2;
    mpreal Rsi, Rsj;

    int N = s.size();

    for (vector<Star>::iterator si = s.begin(); si != s.end()-1; ++si) {
        for (vector<Star>::iterator sj = si+1; sj != s.end(); ++sj) {

            dx = si->r[0] - sj->r[0];
            dy = si->r[1] - sj->r[1];
            dz = si->r[2] - sj->r[2];
            
            dr2 = dx*dx + dy*dy + dz*dz;
            Rsi = this->rmerge * "6"*si->m*c2; // Least stable circular orbit
            Rsj = this->rmerge * "6"*sj->m*c2; // Least stable circular orbit
            
            if (dr2 < std::max(Rsi, Rsj)*std::max(Rsi, Rsj)) {
                cout << "Collision detected!" << endl;
                out.merged = true;
                out.idx0 = distance(begin(s), si);
                out.idx1 = distance(begin(s), sj);
                
                cout << this->merge_idx0 << " " << this->merge_idx1 << endl;
                break;
            }
        }
        if (out.merged) {
            break;
        }
    }
}

vector<double> Cluster::get_data_double() {
    vector<double> ddata;  
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        ddata.push_back(si->m.toDouble());
        ddata.push_back(si->r[0].toDouble());
        ddata.push_back(si->r[1].toDouble());
        ddata.push_back(si->r[2].toDouble());
        ddata.push_back(si->v[0].toDouble());
        ddata.push_back(si->v[1].toDouble());
        ddata.push_back(si->v[2].toDouble());	  
    }	  
    return ddata;
}
vector<mpreal> Cluster::get_data() {
    vector<mpreal> ddata;  
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        ddata.push_back(si->m);
        ddata.push_back(si->r[0]);
        ddata.push_back(si->r[1]);
        ddata.push_back(si->r[2]);
        ddata.push_back(si->v[0]);
        ddata.push_back(si->v[1]);
        ddata.push_back(si->v[2]);	  
    }	  
    return ddata;
}

