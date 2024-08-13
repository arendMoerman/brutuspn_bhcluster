#include <iostream>
using namespace std;

#include "mpCommunicator.h"

int main(int argc, char* argv[]) {

  mpCommunicator comm;
  comm.start_mpi(argc, argv);

  int numBits = 72;
  mpreal::set_default_prec(numBits);
  int numDigits = abs( log10( pow("2", -numBits) ) ).toLong();
  cout.precision(numDigits);
  comm.set_mpreal(numBits/2);

  // Integers
  int a = 0;
  if(comm.get_rank() == 0) a = 1;

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << a << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }

  comm.bcast(a);

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << a << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }

  // bcast mpreal
  mpreal b = "0";
  if(comm.get_rank() == 0) b = acos("-1");

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << b << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }

  comm.bcast(b);

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << b << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }

  // bcast vector<mpreal>
  mpreal b = "0";
  if(comm.get_rank() == 0) b = acos("-1");

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << b << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }

  comm.bcast(b);

  for(int i=0; i<comm.get_size(); i++) {
    if(comm.get_rank() == i) {
      cout << comm.get_rank() << ": " << b << endl;
    }
    for(int j=0; j<1000; j++) {
      mpreal x = acos("-1.2");
    } 
  }
  

  comm.stop_mpi();

  return 0;
}
