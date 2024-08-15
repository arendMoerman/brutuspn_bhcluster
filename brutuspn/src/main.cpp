#include <iostream>
using namespace std;

#include <fstream>
#include <iomanip>

#include <vector>
#include <cmath>
#include <numeric> 
#include <cstdlib>

#include "mpreal.h"
using namespace mpfr;

#include "Timer.h"
#include "Data_Handler.h"

#include "Diagnostics.h"
#include "Initializer.h"

#include "Brutus.h"

void MyOut(int &N, vector<string> &data, ofstream &out) {
    for(int i=0; i<N; i++) {
        for(int j=0; j<7; j++) {
            out << data[i*7+j] << " ";
        }
        out << N << endl;
    }
}

int main(int argc, char* argv[]) {

  int numBits = atoi(argv[7]);
  mpreal::set_default_prec(numBits);
  int numDigits = (int)abs(log10( pow("2.0", -numBits) )).toLong();

  ///////////////////////////////////////////////////////// 
  
  string file_path = "./brutuspn/t.par";
  string value;
  ifstream inFile;
  mpreal t_scale; // timescale in kyr
  inFile.open(file_path);
  if (inFile) {
      while (inFile >> value) {
          t_scale = value;
      }
  }
  inFile.close();

  string outdir = "./outputs/";

  string file_out = outdir + argv[1];

  mpreal t_begin  = argv[2] / t_scale;
  mpreal t_end    = argv[3] / t_scale;
  mpreal dt       = argv[4] / t_scale;
  mpreal eta      = argv[5]; // This one is just given in nbody units by the shell script

  mpreal tolerance = argv[6];

  int nmax = atoi(argv[8]);

  int N		= atoi(argv[9]);
  string config = argv[10];
  vector<string> par;
  for(int i=11; i<argc; i++) {
    par.push_back(argv[i]);
  }

  /////////////////////////////////////////////////////////

  cout.precision(numDigits+16);

  Initializer initializer;
  vector<mpreal> data0 = initializer.generate(N, config, par);
  vector<mpreal> data;
  for(int i=0; i<data0.size(); i++) {
    //cout << data0[i] << endl;
    data.push_back((mpreal)data0[i]);
  //  cout << data[data.size()-1] << endl;
  }
  N = data.size()/7;

  /////////////////////////////////////////////////////////

  Data_Handler data_handler;

  string file_ddata = file_out + ".diag";
  string file_ldata = file_out + ".log";

  ofstream ldata;
  //ddata.precision(numDigits);
  ldata.precision(numDigits);

  //ddata.open(file_ddata.c_str());
  //if(!ddata) {
  //  cerr << "Can't open " << file_ddata << "!" << endl;
  //  return 0;
  //}

  ldata.open(file_ldata.c_str());
  if(!ldata) {
    cerr << "Can't open " << file_ldata << "!" << endl;
    return 0;
  }

  /////////////////////////////////////////////////////////

  mpreal t = t_begin;
  Brutus brutus(t_begin, data, tolerance, numBits, eta, nmax);

  Timer timer;
  double t_cpu = 0;

  cerr << "Simulation started... " << endl;

  vector<string> sdata = brutus.get_data_string();
  vector<mpreal> mdata = brutus.get_data();

  mpreal t0 = (mpreal)atof(argv[2]);
  //data_handler.print(t0, N, t_cpu, data0);
  //data_handler.print(t0, N, t_cpu, data0);

  /////////////////////////////////////////////////////////
  //cout << par[0] << endl;
  cerr << t << "/" << t_end << endl;
  ofstream MyOutfile(file_out + ".out");
  MyOut(N, sdata, MyOutfile);
  
  ofstream Energyfile;
  Energyfile.open(file_out + ".energies");
  Energyfile << std::setprecision(numDigits) << t << " " << brutus.get_energy() << endl;
  
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
  while(t < t_end) {
    timer.start();

    t += dt;
    if(t > t_end) t = t_end;

    bool merge = brutus.evolve(t);
    if(merge) {N--;}

    timer.stop();
    t_cpu += timer.get();

    sdata = brutus.get_data_string();
    mdata = brutus.get_data();
    //data_handler.print(t, N, t_cpu, mdata);
    MyOut(N, sdata, MyOutfile);
    Energyfile << std::setprecision(numDigits) << t << " " << brutus.get_energy() << endl;
    
    ///////////////// analysis stuff ////////////////
    //Ener = brutus.get_Ener();
    //Enerfile << std::setprecision(numDigits) << t;
    //for(vector<mpreal>::iterator e = Ener.begin(); e!=Ener.end(); ++e) {
    //    Enerfile << " " << *e;
    //}
    //Enerfile << endl;
    /////////////////////////////////////////////////
    
    cerr << t << "/" << t_end << endl;    

    //if(sdata[1] == "nan") {
    //    break;
    //}

  }
  Energyfile.close();
  MyOutfile.close();
  //Enerfile.close();

  sdata = brutus.get_data_string();
  mdata = brutus.get_data();
  //data_handler.print(t, N, t_cpu, mdata);

  /////////////////////////////////////////////////////////

  ldata << "N           = " << N << endl;
  ldata << "config      = " << config << endl;
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
  ldata << "x0          = " << sdata[1] << endl;

  cerr << "N           = " << N << endl;
  cerr << "config      = " << config << endl;
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
  cerr << "x0          = " << sdata[1] << endl;

  //ddata.close();
  ldata.close();
  cerr << "Simulation finished!" << endl;

  return 0;
}



    

