// acc, dt, 2-body encounters detect both in double loop

int numNN = 2;
double fiso = 4;

vector<double> data;
int N = data.size()/7;

vector<double> acc(3*N, 0);

double dt = 1e100;

vector< vector<double> > rnn(N)
for(int i=0; i<N; i++) {
  rnn[i].assign(numNN, 1e100);
}

for(int i=0; i<N; i++) {
  for(int j=0; j<N; j++) {
    if(i != j) {
      // Accelerations
      double dx  = data[j*7+1]-data[i*7+1];
      double dy  = data[j*7+2]-data[i*7+2];
      double dz  = data[j*7+3]-data[i*7+3];

      double dr2 = dx*dx + dy*dy + dz*dz;
      double dr = sqrt(dr2);
      double dr3 = dr2*dr;

      double ax = data[j*7]/dr3*dx;
      double ay = data[j*7]/dr3*dy;
      double az = data[j*7]/dr3*dz;

      acc[i*3+0] += ax;
      acc[i*3+1] += ay;
      acc[i*3+2] += az;

      // Timestep size
      double a2 = ax*ax + ay*ay + az*az;
      double dt4 = dr2/a2;
      double dvx = data[j*7+4]-data[i*7+4];
      double dvy = data[j*7+5]-data[i*7+5];
      double dvz = data[j*7+6]-data[i*7+6];
      double dv2 = dvx*dvx + dvy*dvy + dvz*dvz;
      if(dv2 > 0) {
        double dt4_rv = (dr2*dr2)/(dv2*dv2);
        if(dt4_rv < dt4) dt4 = dt4_rv;
      }
      if(dt4 < dt) dt = dt4;

      // Multiples
      for(int k=0; k<numNN; k++) {
        if(dr < rnn[i][k]) {
          for(int q=numNN-1; q>k; q--) {
            rnn[i][k] = rnn[i][k-1];
          }
          rnn[i][k] = rnn;
        }
      }

      ;
    }
  }
}

dt = pow(dt, 0.25);

for(int i=0; i<N; i++) {
  double f = rnn[i][1]/rnn[i][0];
  if(f > fiso) {
    ;
  }
}


