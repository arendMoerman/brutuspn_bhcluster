#include <iostream>
using namespace std;

#include <cstdlib>

#include "Data_Handler.h"
#include "Initializer.h"

int main(int argc, char* argv[]) {
  int N = atoi(argv[1]);
  string config = argv[2];
  vector<string> par;
  for(int i=3; i<argc; i++) {
    par.push_back(argv[i]);
  }

  Data_Handler dh;
  Initializer init;

  vector<double> data = init.generate(N, config, par);
  dh.print(data);
  
  return 0;
}
