#include "Acceleration.h"

// Function to read speed of light from "c.par". Called in constructor
void Acceleration::set_cl() {
    string file_path = "./brutuspn/c.par";
    string value;
    ifstream inFile;
    mpreal c;
    inFile.open(file_path);
    if (inFile) {
        while (inFile >> value) {
            c = value;
        }
    }
    inFile.close();
    this->cl = c;
    this->c2 = "1"/(cl*cl);
    this->c5 = "1"/(cl*cl*cl*cl*cl);
}

// Function to read desired PN terms from "PN.order". Also called in constructor
void Acceleration::set_PN() {
    string file_path = "./brutuspn/PN.order";
    string value;
    ifstream inFile;
    vector<string> PN;
    inFile.open(file_path);
    if (inFile) {
        while (inFile >> value) {
            PN.push_back(value);
        }
    }
    inFile.close();
    PN1p    = stoi(PN[0]);
    PN1c    = stoi(PN[1]);
    PN2     = stoi(PN[2]);
    PN2_5   = stoi(PN[3]);
    PN3     = stoi(PN[4]);
    PN3_5   = stoi(PN[5]);
}

// Setting some handy parameters for the PN terms. Used also in Cluster::energies()
void Acceleration::set_PN_prod(array<mpreal, 3> &vj, array<mpreal, 3> &vi, array<mpreal, 3> &dv, array<mpreal, 3> &dr) {
    vj2     = "0";
    viInvj  = "0";
    dvij2   = "0";
    vijvj   = "0";
    
    vnj     = "0";
    vni     = "0";
    vijn    = "0";
    
    for(int k=0; k<3; k++) {
        vj2     += vj[k]*vj[k];
        viInvj  += vi[k]*vj[k];
        dvij2   += dv[k]*dv[k];
        vijvj   += dv[k]*vj[k];
    
        vnj     += vj[k]*dr[k]/Rij;
        vni     += vi[k]*dr[k]/Rij;
        vijn    += dv[k]*dr[k]/Rij;
    }
}

// 1PN pair
void Acceleration::Acceleration_PN1_pair(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal PairInProduct, Pairn;

    PairInProduct = c2*mj/(Rij*Rij)*("4"*vni - "3"*vnj);
        
    Pairn = c2*mj*apreij*(("4"/Rij*mj + "5"/Rij*mi) - vi2 + 
            "4"*viInvj - "2"*vj2 + "3"*vnj*vnj/"2");
    
    for(int k=0; k<3; k++) {
        da[k] += PairInProduct*dv[k] + Pairn*dr[k];
    }
}

// 1PN cross
void Acceleration::Acceleration_PN1_cross(const Star &si, const Star &sj, array<mpreal, 3> &da, array<mpreal, 3> &dr) {
    array<mpreal, 3> jk;
    array<mpreal, 3> ik;
    array<mpreal, 3> Cross2;
        
    mpreal aprejk, jk2, ik2, Rjk, Rik;
    mpreal mk;
        
    mpreal Cross1i;
    mpreal InProductxx;
    
    for (vector<Star>::iterator sk = s.begin(); sk != s.end(); ++sk) {
        if(*sk == si || *sk == sj) continue;
        
        mk = sk->m;
        // Set InProductxx, jk2 en ik2 to 0 at the beginning of each iteration in sk
        
        jk2         = "0";
        ik2         = "0";
        InProductxx = "0";
        
        for(int k=0; k<3; k++) {
            jk[k] = sj.r[k] - sk->r[k];
            ik[k] = si.r[k] - sk->r[k];
            jk2 += jk[k]*jk[k];
            ik2 += ik[k]*ik[k];
            InProductxx += dr[k]*jk[k];
        }
        
        Rjk = sqrt(jk2 + eps2);
        Rik = sqrt(ik2 + eps2);
        
        aprejk = "1"/(Rjk*Rjk*Rjk);

        Cross1i += mk*(("1"/Rjk + "4"/Rik) - "1"*aprejk*InProductxx/"2"); 
        
        for(int k=0; k<3; k++) {
            Cross2[k] += mk*aprejk*jk[k];
        }
    }
    
    for(int k=0; k<3; k++) {
        da[k] += c2*mj*(apreij*dr[k]*Cross1i - "3.5"/Rij*Cross2[k]);
    }
}

// 2PN
void Acceleration::Acceleration_PN2(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor2n, factor2v;
    
    factor2n = c2*c2*apreij*mj*("-2"*vj2*vj2 + ("4"*vj2 - "2"*viInvj)*viInvj + ("1.5"*vi2 + "4.5"*vj2 - "6"*viInvj - 
            "1.875"*vnj*vnj)*vnj*vnj - ("14.25"*mi*mi + "9"*mj*mj + "34.5"*mi*mj)/(Rij*Rij) + 
        mi/Rij*("-3.75"*vi2 + "1.25"*vj2 - "2.5"*viInvj + "19.5"*vni*vni - "39"*vni*vnj + "8.5"*vnj*vnj) + 
        mj/Rij*("4"*vj2 - "8"*viInvj + "2"*vni*vni - "4"*vni*vnj - "6"*vnj*vnj));
      
    factor2v = c2*c2*apreij*Rij*mj*(mi/Rij*("-15.75"*vni + "13.75"*vnj) - "2"*mj/Rij*(vni + vnj) + vi2*vnj + "4"*vj2*vni -
            "5"*vj2*vnj - "4"*viInvj*vijn - "6"*vni*vnj*vnj + "4.5"*vnj*vnj*vnj);
    
    for(int k=0; k<3; k++) {
        da[k] += factor2n*dr[k] + factor2v*dv[k];
    }
}

// 2.5PN
void Acceleration::Acceleration_PN2_5(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor2_5n, factor2_5v;
    
    factor2_5n = "0.8"*c5*mi*mj*apreij*vijn/Rij*("-6"*mi/Rij + "52"*mj/(Rij*"3") + "3"*dvij2);
    
    factor2_5v = "0.8"*c5*mi*mj*apreij*("2"*mi/Rij - "8"*mj/Rij - dvij2);
    
    for(int k=0; k<3; k++) {
        da[k] += factor2_5n*dr[k] + factor2_5v*dv[k];
    }
}

// 3PN
void Acceleration::Acceleration_PN3(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor3n, factor3v;
    
    factor3n = c2*c2*c2*apreij*mj*(("2.1875"*vnj*vnj*vnj*vnj - "1.875"*vnj*vnj*vi2 + "7.5"*vnj*vnj*viInvj + "3"*viInvj*viInvj - 
            "7.5"*vnj*vnj*vj2 + "1.5"*vi2*vj2 - "12"*viInvj*vj2 + "7.5"*vj2*vj2)*vnj*vnj - "2"*viInvj*viInvj*vj2 + "4"*viInvj*vj2*vj2 - 
            "2"*vj2*vj2*vj2 +
        mi/Rij*("-21.375"*vni*vni*vni*vni + "85.5"*vni*vni*vni*vnj - "180.75"*vni*vni*vnj*vnj + "191.5"*vni*vnj*vnj*vnj -
            "56.875"*vnj*vnj*vnj*vnj + "57.25"*vni*vni*vi2 - "102.5"*vni*vnj*vi2 + "47.75"*vnj*vnj*vi2 - "11.375"*vi2*vi2 - "114.5"*vni*vni*viInvj +
            "244"*vni*vnj*viInvj - "127.5"*vnj*vnj*viInvj + "45.5"*vi2*viInvj - "44.25"*viInvj*viInvj + "57.25"*vni*vni*vj2 - "141.5"*vni*vnj*vj2 +
            "64.75"*vnj*vnj*vj2 - "22.75"*vi2*vj2 + "43"*viInvj*vj2 - "10.125"*vj2*vj2) + 
        mj/Rij*("-6"*vni*vni*vnj*vnj + "12"*vni*vnj*vnj*vnj + "6"*vnj*vnj*vnj*vnj + "4"*vni*vnj*viInvj + "12"*vnj*vnj*viInvj +
            "4"*viInvj*viInvj - "4"*vni*vnj*vj2 - "12"*vnj*vnj*vj2 - "8"*viInvj*vj2 + "4"*vj2*vj2) +
        mj*mj/(Rij*Rij)*("2"*vni*vnj - vni*vni + "21.5"*vnj*vnj + "18"*viInvj - "9"*vj2) +
        mi*mj/(Rij*Rij)*("51.875"*vni*vni - "93.75"*vni*vnj + "139.125"*vnj*vnj - "9.609375"*pi*pi*vijn*vijn + "18"*vi2 + "1.921875"*pi*pi*dvij2 +
            "33"*viInvj - "16.5"*vj2) +
        mi*mi/(Rij*Rij)*("-258.625"*vni*vni + "543"*vni*vnj - "234.75"*vnj*vnj + "58.875"*vi2 - "89.25"*viInvj + "44.625"*vj2) +
            ("16"*mj*mj*mj + mi*mi*mj*("547" - "123"*pi*pi/"16")/"3" - "13"*mi*mi*mi/"12" + mi*mj*mj*("545" - "123"*pi*pi/"16")/"3")*apreij);
        
    factor3v = c2*c2*c2*mj/(Rij*Rij)*("7.5"*vni*vnj*vnj*vnj*vnj - "5.625"*vnj*vnj*vnj*vnj*vnj - "1.5"*vnj*vnj*vnj*vi2 +
            "6"*vni*vnj*vnj*viInvj - "6"*vnj*vnj*vnj*viInvj - "2"*vnj*viInvj*viInvj - "12"*vni*vnj*vnj*vj2 + "12"*vnj*vnj*vnj*vj2 +
            vnj*vi2*vj2 - "4"*vni*viInvj*vj2 + "8"*vnj*viInvj*vj2 + "4"*vni*vj2*vj2 - "7"*vnj*vj2*vj2 +
        mj/Rij*("-2"*vni*vni*vnj + "8"*vni*vnj*vnj + "2"*vnj*vnj*vnj + "2"*vni*viInvj + "4"*vnj*viInvj - "2"*vni*vj2 - "4"*vnj*vj2) +
        mi/Rij*("-60.75"*vni*vni*vni + "141.25"*vni*vni*vnj - "67.25"*vni*vnj*vnj - "95"*vnj*vnj*vnj/"12" + "25.875"*vni*vi2 -
            "17.125"*vnj*vi2 - "36"*vni*viInvj + "6.75"*vnj*viInvj + "10.125"*vni*vj2 + "10.375"*vnj*vj2) +
        (mj*mj*("4"*vni + "5"*vnj) +
        mi*mj*("-38.375"*vni + "59.875"*vnj + "3.84375"*pi*pi*vijn) +
        mi*mi*("77.75"*vni - "89.25"*vnj))/(Rij*Rij));
    
    for(int k=0; k<3; k++) {
        da[k] += factor3n*dr[k] + factor3v*dv[k];
    }
}

// 3.5PN
void Acceleration::Acceleration_PN3_5(array<mpreal, 3> &da, array<mpreal, 3> &dr, array<mpreal, 3> &dv) {
    mpreal factor3_5n, factor3_5v;
    
    factor3_5n = c2*c5*apreij/Rij*(mi*mi*mi*mj*apreij*("3992"*vni/"105" - "4328"*vnj/"105") +
        apreij*mi*mi*mj*mj*("-13576"*vni/"105" + "2872"*vnj/"21") - "3172"*apreij*mi*mj*mj*mj*vijn/"21" +
        mi*mi*mj/Rij*("48"*vni*vni*vni - "139.2"*vni*vni*vnj + "148.8"*vni*vnj*vnj - "57.6"*vnj*vnj*vnj - 
            "4888"*vni*vi2/"105" + "5056"*vnj*vi2/"105" + "2056"*vni*viInvj/"21" - "2224"*vnj*viInvj/"21" -
            "1028"*vni*vj2/"21" + "5812"*vnj*vj2/"105") +
        mi*mj*mj/Rij*("-116.4"*vni*vni*vni + "349.2"*vni*vni*vnj - "390.8"*vni*vnj*vnj +
            "158"*vnj*vnj*vnj + "3568"*vijn*vi2/"105" - "2864"*vni*viInvj/"35" + "10048"*vnj*viInvj/"105" +
            "1432"*vni*vj2/"35" - "5752"*vnj*vj2/"105") +
        mi*mj*("-56"*vijn*vijn*vijn*vijn*vijn + "60"*vni*vni*vni*dvij2 - "180"*vni*vni*vnj*dvij2 +
            "174"*vni*vnj*vnj*dvij2 - "54"*vnj*vnj*vnj*dvij2 - "246"*vijn*vi2*vi2/"35" +
            "1068"*vni*vi2*viInvj/"35" - "984"*vnj*vi2*viInvj/"35" - "1068"*vni*viInvj*viInvj/"35" +
            "180"*vnj*viInvj*viInvj/"7" - "534"*vni*vi2*vj2/"35" + "90"*vnj*vi2*vj2/"7" + "984"*vni*viInvj*vj2/"35" -
            "732"*vnj*viInvj*vj2/"35" - "204"*vni*vj2*vj2/"35" + "24"*vnj*vj2*vj2/"7"));
             
    factor3_5v = c2*c5*apreij*("-184"*mi*mi*mi*mj/("21"*Rij*Rij) + "6224"*mi*mi*mj*mj*apreij/"105" +
            "6388"*apreij*mi*mj*mj*mj/"105" +
        mi*mi*mj/Rij*("52"*vni*vni/"15" - "56"*vni*vnj/"15" - "44"*vnj*vnj/"15" - "132"*vi2/"35" + "152"*viInvj/"35" -
            "48"*vj2/"35") +
        mi*mj*mj/Rij*("454"*vni*vni/"15" - "74.4"*vnj*vni + "854"*vnj*vnj/"15" - "152"*vi2/"21" + 
            "2864"*viInvj/"105" - "1768"*vj2/"105") +
        mi*mj*("60"*vijn*vijn*vijn*vijn - "69.6"*vni*vni*dvij2 + "136.8"*vni*vnj*dvij2 -
            "66"*vnj*vnj*dvij2 + "334"*vi2*vi2/"35" - "1336"*vi2*viInvj/"35" + "1308"*viInvj*viInvj/"35" +
            "654"*vi2*vj2/"35" - "1252"*viInvj*vj2/"35" + "292"*vj2*vj2/"35"));
    
    for(int k=0; k<3; k++) {
        da[k] += factor3_5n*dr[k] + factor3_5v*dv[k];
    }
}

// Function to calc barycenter
void Acceleration::get_bary() {
    
    mpreal M = "0";
    barv.fill("0");
    bara.fill("0");
    for(vector<Star>::iterator si = s.begin(); si != s.end(); ++si)  {
        for(int k=0; k<3; k++) {

            M += si->m;
            barv[k] += si->m*si->v[k];
            bara[k] += si->m*si->a[k];
        }
    }
    for(int k=0; k<3; k++) {
        barv[k] /= M;
        bara[k] /= M;
    }
    //cout << barv[0] << " " << barv[1] << endl;
}

// Here the acceleration used during Bulirsch-Stoer steps
void Acceleration::calcAcceleration() {
    array<mpreal, 3> vi;
    array<mpreal, 3> dr;
    array<mpreal, 3> dv;
    array<mpreal, 3> da;
    mpreal dr2 = "0";

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->a.fill("0");
    }

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        array<mpreal, 3> vj;
        vi2 = "0";
        
        for(int k=0; k<3; k++) {
            vi[k] = si->v[k];
            vi2  += vi[k]*vi[k];
        }
        mi = si->m;
        
        for (vector<Star>::iterator sj = s.begin(); sj != s.end(); ++sj) {
            if(*sj == *si) continue;
            
            dr2 = "0";
            da.fill("0");
            
            for(int k=0; k<3; k++) {
                vj[k] = sj->v[k];
                dr[k] = si->r[k]-sj->r[k];
                dv[k] = si->v[k]-sj->v[k];
                dr2  += dr[k]*dr[k];
            }
            mj = sj->m;

            Rij = sqrt(dr2 + eps2);
            apreij = "1"/(Rij*Rij*Rij);

            set_PN_prod(vj, vi, dv, dr);

            if (PN1p)   Acceleration_PN1_pair(da, dr, dv);
            if (PN1c)   Acceleration_PN1_cross(*si, *sj, da, dr);
            if (PN2)    Acceleration_PN2(da, dr, dv);
            if (PN2_5)  Acceleration_PN2_5(da, dr, dv);
            if (PN3)    Acceleration_PN3(da, dr, dv);
            if (PN3_5)  Acceleration_PN3_5(da, dr, dv);

            for(int k=0; k<3; k++) {
                da[k] -= mj*apreij*dr[k];
                si->a[k] += da[k];
            }
        }
    }
}

// Here, the acceleration calculated at the very beginning of a new time step.
// The initial stepsize is calculated here
void Acceleration::calcAcceleration_dt() {
    array<mpreal, 3> vi;
    array<mpreal, 3> dr;
    array<mpreal, 3> dv;
    array<mpreal, 3> da;
    mpreal da2 = "0";
    mpreal dr2 = "0";

    dt = "1e100";
    mpreal mydt = "0";
    
    int N = s.size();
    
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        si->a.fill("0");
    }

    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        array<mpreal, 3> vj;
        vi2 = "0";
        
        for(int k=0; k<3; k++) {
            vi[k] = si->v[k];
            vi2  += vi[k]*vi[k];
        }
        mi = si->m;
        
        for (vector<Star>::iterator sj = s.begin(); sj != s.end(); ++sj) {
            if(*sj == *si) continue;
            
            dr2 = "0";
            da.fill("0");
            for(int k=0; k<3; k++) {
                vj[k] = sj->v[k];
                dr[k] = si->r[k]-sj->r[k];
                dv[k] = si->v[k]-sj->v[k];
                dr2  += dr[k]*dr[k];
            }
            mj = sj->m;

            Rij = sqrt(dr2 + eps2);
            apreij = "1"/(Rij*Rij*Rij);

            set_PN_prod(vj, vi, dv, dr);

            if (PN1p)   Acceleration_PN1_pair(da, dr, dv);
            if (PN1c)   Acceleration_PN1_cross(*si, *sj, da, dr);
            if (PN2)    Acceleration_PN2(da, dr, dv);
            if (PN2_5)  Acceleration_PN2_5(da, dr, dv);
            if (PN3)    Acceleration_PN3(da, dr, dv);
            if (PN3_5)  Acceleration_PN3_5(da, dr, dv);
            
            da2 = "0";
            for(int k=0; k<3; k++) {
                da[k]    -= mj*apreij*dr[k];
                si->a[k] += da[k];
                da2      += da[k]*da[k];
            }
            mydt = dr2 / da2;
            if(mydt < dt) dt = mydt;
        }
    }
    dt = sqrt(sqrt(dt)); //pow(dt, "0.25");
}

// Toy routine for getting acceleration output
vector<array<mpreal, 3>> Acceleration::getAcceleration() {
    int N = s.size();
    
    vector<array<mpreal, 3>> ai;
    for (vector<Star>::iterator si = s.begin(); si != s.end(); ++si) {
        ai.push_back(si->a);
    }
    return ai;
}
