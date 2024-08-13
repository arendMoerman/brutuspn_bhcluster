#include "Tools.h"

vector< vector<int> > Tools::get_nn(vector<double> &data, int numNN) {
  int N = data.size()/7;

  vector< vector<int> > nn_index(N);
  vector< vector<double> > nn_dr2(N);
  for(int i=0; i<N; i++) {
    nn_index[i].assign(numNN, 0);
    nn_dr2[i].assign(numNN, 1e100);
  }

  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      double dx = data[j*7+1]-data[i*7+1];
      double dy = data[j*7+2]-data[i*7+2];
      double dz = data[j*7+3]-data[i*7+3];
      double dr2 = dx*dx + dy*dy + dz*dz;

      for(int k=0; k<numNN; k++) {
        if(dr2 < nn_dr2[i][k]) {
          for(int q=numNN-1; q>k; q--) {
	    nn_dr2[i][q] = nn_dr2[i][q-1];
	    nn_index[i][q] = nn_index[i][q-1];
          }
          nn_dr2[i][k] = dr2;
          nn_index[i][k] = j;
          break;
        }  
      }

      for(int k=0; k<numNN; k++) {
        if(dr2 < nn_dr2[j][k]) {
          for(int q=numNN-1; q>k; q--) {
	    nn_dr2[j][q] = nn_dr2[j][q-1];
	    nn_index[j][q] = nn_index[j][q-1];
          }
          nn_dr2[j][k] = dr2;
          nn_index[j][k] = i;
          break;
        }  
      }
      
    }
  }

  return nn_index;
}  







