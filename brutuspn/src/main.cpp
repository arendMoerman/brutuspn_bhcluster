#include <iostream>
using namespace std;

#include <fstream>
#include <iomanip>

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>
#include <array>

#include "mpreal.h"
using namespace mpfr;

#include "Timer.h"
#include "Data_Handler.h"

#include "Diagnostics.h"
#include "Initializer.h"

#include "Brutus.h"

void MyOut(int &N, vector<string> &data, ofstream &out, vector<mpreal> &mrvt) {
    mpreal factor;
    for(int i=0; i<N; i++) {
        factor = mrvt[0]; // Write mass
        for(int j=0; j<13; j++) {
            if(j > 0) {factor = mrvt[1];} // Write positions
            if(j > 3) {factor = mrvt[2];} // Write velocities
            if(j > 6) {factor = mrvt[2]/mrvt[4];} // Write accelerations
            if(j > 9) {factor = mrvt[2]/mrvt[4]/mrvt[4];} // Write jerk
            out << data[i*13+j] * factor << " ";
        }
        out << N << endl;
    }
}

int main(int argc, char* argv[]) {

  int numBits = atoi(argv[7]);
  mpreal::set_default_prec(numBits);
  int numDigits = (int)abs(log10( pow("2.0", -numBits) )).toLong();

  ///////////////////////////////////////////////////////// 

  vector<mpreal> mrvt;

  string file_path = "./brutuspn/mrvt.par";
  string value;
  ifstream inFile;
  inFile.open(file_path);
  if (inFile) {
      while (inFile >> value) {
          mrvt.push_back(value);
      }
  }
  inFile.close();

  // Also add time scaling in seconds at end of mrvt vector, for acc and jerk scaling
  mrvt.push_back(mrvt[3] * 1000 * 365 * 24 * 60 * 60);

  string outdir = "./outputs/";

  string file_out = outdir + argv[1];

  mpreal t_begin  = argv[2] / mrvt[3];
  mpreal t_end    = argv[3] / mrvt[3];
  mpreal dt       = argv[4] / mrvt[3];
  mpreal eta      = argv[5]; // This one is just given in nbody units by the shell script

  int tol_power = atoi(argv[6]);

  mpreal tolerance = pow(10, -1*tol_power);

  int nmax = atoi(argv[8]);

  int N		= atoi(argv[9]);

  mpreal r_merge = argv[10];

  string config = argv[11];
  vector<string> par;
  for(int i=12; i<argc; i++) {
    par.push_back(argv[i]);
  }

  /////////////////////////////////////////////////////////

  cout.precision(numDigits+16);

  Initializer initializer;
  vector<mpreal> data0 = initializer.generate(N, config, par);
  vector<mpreal> data;
  for(int i=0; i<data0.size(); i++) {
    data.push_back((mpreal)data0[i]);
  }
  N = data.size()/7;

  /////////////////////////////////////////////////////////

  Data_Handler data_handler;

  string file_ddata = file_out + ".diag";
  string file_ldata = file_out + ".log";

  ofstream ldata;
  ldata.precision(numDigits);

  ldata.open(file_ldata.c_str());
  if(!ldata) {
    cerr << "Can't open " << file_ldata << "!" << endl;
    return 0;
  }

  /////////////////////////////////////////////////////////

  mpreal t = t_begin;
  Brutus brutus(t_begin, data, tolerance, numBits, eta, nmax, r_merge);

  Timer timer;
  double t_cpu = 0;

  cerr << "Simulation started... " << endl;

  vector<string> sdata = brutus.get_data_string();
  //vector<mpreal> mdata = brutus.get_data();

  mpreal t0 = (mpreal)atof(argv[2]);
  //data_handler.print(t0, N, t_cpu, data0);
  //data_handler.print(t0, N, t_cpu, data0);

  /////////////////////////////////////////////////////////
  //cout << par[0] << endl;
  //cerr << t << "/" << t_end << endl;
  ofstream MyOutfile(file_out + ".out");
  MyOut(N, sdata, MyOutfile, mrvt);
  
  ofstream Energyfile;
  Energyfile.open(file_out + ".energies");
  Energyfile << std::setprecision(numDigits) << t*mrvt[3] << " " << brutus.get_energy() << endl;
  
  //ofstream Enerfile;
  //Enerfile.open("Ener.out");
  //vector<mpreal> Ener;
  /*
  vector<mpreal> Ener = brutus.get_Ener();
  Enerfile << std::setprecision(numDigits) << t;
  for(vector<mpreal>::iterator e = Ener.begin(); e!=Ener.end(); ++e) {
    Enerfile << " " << *e;
  }
  Enerfile << endl;
  */
  float prog;
  mergerOut merge;

  while(t < t_end) {
    timer.start();

    t += dt;
    if(t > t_end) t = t_end;

    brutus.evolve(t, merge);
    if(merge.merged) {
        break;
    }

    timer.stop();
    t_cpu += timer.get();

    sdata = brutus.get_data_string();
    //mdata = brutus.get_data();
    MyOut(N, sdata, MyOutfile, mrvt);
    Energyfile << std::setprecision(numDigits) << t*mrvt[3] << " " << brutus.get_energy() << endl;
    
    ///////////////// analysis stuff ////////////////
    //Ener = brutus.get_Ener();
    //Enerfile << std::setprecision(numDigits) << t;
    //for(vector<mpreal>::iterator e = Ener.begin(); e!=Ener.end(); ++e) {
    //    Enerfile << " " << *e;
    //}
    //Enerfile << endl;
    /////////////////////////////////////////////////

    prog = (float)(t / t_end) * 100;

    printf("Progress: %f / 100\r", prog);
    fflush(stdout);
    //cerr << t << "/" << t_end << endl;    

    //if(sdata[1] == "nan") {
    //    break;
    //}

  }
  Energyfile.close();
  MyOutfile.close();
  //Enerfile.close();

  sdata = brutus.get_data_string();
  //mdata = brutus.get_data();
  //data_handler.print(t, N, t_cpu, mdata);

  /////////////////////////////////////////////////////////

  ldata << "N           = " << N << endl;
  ldata << "param       = ";
  for(int i=0; i<par.size(); i++) ldata << par[i] << " ";
  ldata << endl;
  ldata << "t_begin     = " << t_begin << endl;
  ldata << "t_end       = " << t << endl; 
  ldata << "dt          = " << dt << endl;
  ldata << "eta         = " << eta << endl;
  ldata << "e           = " << tolerance << endl;
  ldata << "Lw          = " << numBits << endl;
  ldata << "t_cpu       = " << t_cpu << endl;
  ldata << "mergers     = " << merge.idx0 << " " << merge.idx1 << endl;

  cerr << "N           = " << N << endl;
  cerr << "param       = ";
  for(int i=0; i<par.size(); i++) cerr << par[i] << " ";
  cerr << endl;
  cerr << "t_begin     = " << t_begin << endl;
  cerr << "t_end       = " << t << endl; 
  cerr << "dt          = " << dt << endl;
  cerr << "eta         = " << eta << endl;
  cerr << "e           = " << tolerance << endl;
  cerr << "Lw          = " << numBits << endl;
  cerr << "t_cpu       = " << t_cpu << endl;
  cerr << "mergers     = " << merge.idx0 << " " << merge.idx1 << endl;

  //ddata.close();
  ldata.close();
  cerr << "Simulation finished!" << endl;

  return 0;
}

