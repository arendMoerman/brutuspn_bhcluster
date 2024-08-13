#include "mpCommunicator.h"

//////////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////////
mpCommunicator::mpCommunicator() {
  ;
}
//////////////////////////////////////////////////////////////////////////
// MPI
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::start_mpi(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  name = processor_name;  
  new_comm = MPI_COMM_WORLD;
}
void mpCommunicator::stop_mpi() {
  MPI_Finalize();
}

void mpCommunicator::make_two_groups() {
  /* Extract the original group handle */
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

  /* Divide tasks into two distinct groups based upon rank */
  vector<int> v1, v2;
  for(int i=0; i<size/2; i++) v1.push_back(i);
  for(int i=size/2; i<size; i++) v2.push_back(i);
  int* ranks1 = &v1[0];
  int* ranks2 = &v2[0];
  MPI_Group_incl(orig_group, size/2, ranks1, &new_group1);
  MPI_Group_incl(orig_group, size/2, ranks2, &new_group2);

  /* Create new communicators */
  MPI_Comm_create(MPI_COMM_WORLD, new_group1, &new_comm1);
  MPI_Comm_create(MPI_COMM_WORLD, new_group2, &new_comm2);

  /* Keep communicators different than MPI_COMM_NULL based on the rank */
  if(rank < size/2) {
    new_comm = new_comm1;
    group = 1;
  }
  else {
    new_comm = new_comm2;
    group = 2;
  }

  /* Get new ranks */
  MPI_Comm_rank(new_comm, &rank);
  MPI_Comm_size(new_comm, &size);
}
void mpCommunicator::free_groups() {
  /* Free grouping */
  MPI_Comm_free(&new_comm);
  MPI_Group_free(&new_group1);
  MPI_Group_free(&new_group2);
}
//////////////////////////////////////////////////////////////////////////
// Work
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::divide_work(int numElements) {
  numStar = numElements;
  int myNumber = (int)(numStar/size);
  int restNumber = numStar - myNumber*size;
  for(int i=0; i<restNumber; i++) {
    if(i == rank) myNumber++;
  }
  gather(myNumber, all_numbers, new_comm);
  bcast(all_numbers, new_comm);
  my_begin = 0;
  for(int i=0; i<rank; i++) {
    my_begin += all_numbers[i];
  }
  my_end = my_begin + myNumber;
  gather(my_begin, all_begins, new_comm);
  bcast(all_begins, new_comm);
}

//////////////////////////////////////////////////////////////////////////
// Mpreal
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::set_mpreal(int numDigits) {
  this->numDigits = numDigits;
}
void mpCommunicator::format(string &x) {
  int flag = 0;
  string::iterator it;
  while(flag == 0 && x.size() > 1) {
    flag = 1;
    if(x[0] == '0') {
      it = x.begin();
      x.erase (it);
      flag = 0;
    }
  }
}
void mpCommunicator::convert(vector<char> &x, mpreal &y) {
  x.clear();
  x.resize(numDigits);
  string z = y.toString();
  if(z.size() < numDigits) {
    string addon;
    for(int i=0; i<numDigits-z.size(); i++) {
      addon.push_back('0');
    }
    z = addon + z;
  }
  for(int i=0; i<numDigits; i++) {
    x[i] = z[i];
  }
}
void mpCommunicator::deconvert(mpreal &x, vector<char> &y) {
  string z;
  for(int i=0; i<numDigits; i++) {
    z.push_back( y[i] );
  }
  format(z);
  x = z.c_str();
}
void mpCommunicator::serialize(vector<char> &y, vector<mpreal> &x) {
  int N = x.size();
  int Q = N*numDigits;
  y.resize(Q);
  vector<string> z(N);
  for(int i=0; i<N; i++) {
    z[i] = x[i].toString();
  }
  for(int i=0; i<N; i++) {
    if(z[i].size() < numDigits) {
      string addon;
      for(int j=0; j<numDigits-z[i].size(); j++) {
        addon.push_back('0');
      }
      z[i] = addon + z[i];
    }
    else if(z[i].size() > numDigits) {
      size_t found;
      string str = "e";	
      found = z[i].find(str);
      if(found != string::npos) {
	int numErase = z[i].size()-numDigits;
	z[i].erase (int(found)-numErase, numErase);
      }	
    }
  }  
  int k=0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<numDigits; j++) {
      y[k] = z[i][j];
      k++;
    }
  }
}
void mpCommunicator::deserialize(vector<mpreal> &x, vector<char> &y) {
  int N = y.size()/numDigits;
  x.clear();
  x.resize(N);
  vector<string> z(N);
  int k = 0;
  for(int i=0; i<N; i++) {
    string temp;
    for(int j=0; j<numDigits; j++) {
      temp.push_back( y[k] );
      k++;
    }
    z[i] = temp;
  }
  for(int i=0; i<N; i++) {
    format(z[i]);
    x[i] = z[i].c_str();
  }  
}
//////////////////////////////////////////////////////////////////////////
// Getters
//////////////////////////////////////////////////////////////////////////
string mpCommunicator::get_name() {
  return name;
}

int mpCommunicator::get_rank() {
  return rank;
}
int mpCommunicator::get_size() {
  return size;
}
int mpCommunicator::get_my_begin() {
  return my_begin;
}
int mpCommunicator::get_my_end() {
  return my_end;
}

int mpCommunicator::get_group() {
  return group;
}
MPI_Comm mpCommunicator::get_comm_world() {
  return MPI_COMM_WORLD;
}
MPI_Comm mpCommunicator::get_comm_group1() {
  return new_comm1;
}
MPI_Comm mpCommunicator::get_comm_group2() {
  return new_comm2;
}
//////////////////////////////////////////////////////////////////////////
// INT
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::bcast(int &x) {
  MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);   
}
void mpCommunicator::bcast(vector<int> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_INT, 0, MPI_COMM_WORLD);
}
void mpCommunicator::gather(int &x, vector<int> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_INT, &y.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);  
}
void mpCommunicator::gather(vector<int> &x, vector<int> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_INT, &y.front(), &N_col.front(), &i_col.front(), MPI_INT, 0, MPI_COMM_WORLD); 
}
void mpCommunicator::join(vector<int> &x, vector<int> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_INT, &y.front(), &all_number.front(), &all_begin.front(), MPI_INT, 0, MPI_COMM_WORLD);
}
void mpCommunicator::send_id() {
  int send_msg = this->rank;
  int send_tag = this->rank;
  MPI_Send(&send_msg, 1, MPI_INT, 0, send_tag, MPI_COMM_WORLD);
}
void mpCommunicator::recv_id(int &src_id) {
  MPI_Status status;
  int recv_msg;
  MPI_Recv(&recv_msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  src_id = status.MPI_SOURCE;
}
void mpCommunicator::send_ic(int &dest_id, int &ic_index) {
  MPI_Send(&ic_index, 1, MPI_INT, dest_id, 0, MPI_COMM_WORLD);
}
void mpCommunicator::recv_ic(int &ic_index) {
  MPI_Status status;
  MPI_Recv(&ic_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
}

void mpCommunicator::bcast(int &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_INT, 0, comm);   
}
void mpCommunicator::bcast(vector<int> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_INT, 0, comm);
}
void mpCommunicator::gather(int &x, vector<int> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_INT, &y.front(), 1, MPI_INT, 0, comm);  
}
void mpCommunicator::gather(vector<int> &x, vector<int> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_INT, &y.front(), &N_col.front(), &i_col.front(), MPI_INT, 0, comm); 
}
void mpCommunicator::join(vector<int> &x, vector<int> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_INT, &y.front(), &all_number.front(), &all_begin.front(), MPI_INT, 0, comm);
}
void mpCommunicator::send_id(MPI_Comm comm) {
  int send_msg = this->rank;
  int send_tag = this->rank;
  MPI_Send(&send_msg, 1, MPI_INT, 0, send_tag, comm);
}
void mpCommunicator::recv_id(int &src_id, MPI_Comm comm) {
  MPI_Status status;
  int recv_msg;
  MPI_Recv(&recv_msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
  src_id = status.MPI_SOURCE;
}
void mpCommunicator::send_ic(int &dest_id, int &ic_index, MPI_Comm comm) {
  MPI_Send(&ic_index, 1, MPI_INT, dest_id, 0, comm);
}
void mpCommunicator::recv_ic(int &ic_index, MPI_Comm comm) {
  MPI_Status status;
  MPI_Recv(&ic_index, 1, MPI_INT, 0, 0, comm, &status);
}
/////////////////////////////////////////////////////////////////////////
// FLOAT
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::bcast(float &x) {
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);   
}
void mpCommunicator::bcast(vector<float> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_FLOAT, 0, MPI_COMM_WORLD);
}
void mpCommunicator::gather(float &x, vector<float> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_FLOAT, &y.front(), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);  
}
void mpCommunicator::gather(vector<float> &x, vector<float> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_FLOAT, &y.front(), &N_col.front(), &i_col.front(), MPI_FLOAT, 0, MPI_COMM_WORLD); 
}
void mpCommunicator::join(vector<float> &x, vector<float> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_FLOAT, &y.front(), &all_number.front(), &all_begin.front(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void mpCommunicator::bcast(float &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, comm);   
}
void mpCommunicator::bcast(vector<float> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_FLOAT, 0, comm);
}
void mpCommunicator::gather(float &x, vector<float> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_FLOAT, &y.front(), 1, MPI_FLOAT, 0, comm);  
}
void mpCommunicator::gather(vector<float> &x, vector<float> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_FLOAT, &y.front(), &N_col.front(), &i_col.front(), MPI_FLOAT, 0, comm); 
}
void mpCommunicator::join(vector<float> &x, vector<float> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_FLOAT, &y.front(), &all_number.front(), &all_begin.front(), MPI_FLOAT, 0, comm);
}
/////////////////////////////////////////////////////////////////////////
// DOUBLE
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::bcast(double &x) {
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
}
void mpCommunicator::bcast(vector<double> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void mpCommunicator::gather(double &x, vector<double> &y) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_DOUBLE, &y.front(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
}
void mpCommunicator::gather(vector<double> &x, vector<double> &y) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_DOUBLE, &y.front(), &N_col.front(), &i_col.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
}
void mpCommunicator::join(vector<double> &x, vector<double> &y) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_DOUBLE, &y.front(), &all_number.front(), &all_begin.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void mpCommunicator::bcast(double &x, MPI_Comm comm) {
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, comm);   
}
void mpCommunicator::bcast(vector<double> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  MPI_Bcast(&x.front(), N, MPI_DOUBLE, 0, comm);
}
void mpCommunicator::gather(double &x, vector<double> &y, MPI_Comm comm) {
  y.resize(size);
  MPI_Gather(&x, 1, MPI_DOUBLE, &y.front(), 1, MPI_DOUBLE, 0, comm);  
}
void mpCommunicator::gather(vector<double> &x, vector<double> &y, MPI_Comm comm) {
  int N = x.size();
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  y.resize(N_sum);
  MPI_Gatherv(&x.front(), N, MPI_DOUBLE, &y.front(), &N_col.front(), &i_col.front(), MPI_DOUBLE, 0, comm); 
}
void mpCommunicator::join(vector<double> &x, vector<double> &y, MPI_Comm comm) {
  int M = x.size();
  y.resize(M);
  int N = M/numStar;
  int begin = my_begin*N;
  int number = all_numbers[rank]*N;
  vector<int> all_begin = all_begins;
  vector<int> all_number = all_numbers;
  for(int i=0; i<size; i++) {
    all_begin[i] *= N;
    all_number[i] *= N;
  }
  MPI_Gatherv(&x[begin], number, MPI_DOUBLE, &y.front(), &all_number.front(), &all_begin.front(), MPI_DOUBLE, 0, comm);
}
/////////////////////////////////////////////////////////////////////////
// MPREAL
//////////////////////////////////////////////////////////////////////////
void mpCommunicator::bcast(mpreal &x) {
  vector<char> y(numDigits);
  convert(y, x);
  MPI_Bcast(&y.front(), numDigits, MPI_CHAR, 0, MPI_COMM_WORLD);  
  deconvert(x, y);
}
void mpCommunicator::bcast(vector<mpreal> &x) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  vector<char> y(N*numDigits);
  serialize(y, x);  
  MPI_Bcast(&y.front(), N*numDigits, MPI_CHAR, 0, MPI_COMM_WORLD);   
  deserialize(x, y);
}
void mpCommunicator::gather(mpreal &x, vector<mpreal> &y) {
  y.clear();
  y.resize(size);
  vector<char> z(numDigits);
  convert(z, x);
  vector<char> col(size*numDigits);
  MPI_Gather( &z.front(), numDigits, MPI_CHAR, &col.front(), numDigits, MPI_CHAR, 0, MPI_COMM_WORLD);
  deserialize(y, col); 
}
void mpCommunicator::gather(vector<mpreal> &x, vector<mpreal> &y) {
  int N = x.size()*numDigits;
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  vector<char> xx;
  serialize(xx, x);
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  vector<char> yy(N_sum);     
  MPI_Gatherv(&xx.front(), N, MPI_CHAR, &yy.front(), &N_col.front(), &i_col.front(), MPI_CHAR, 0, MPI_COMM_WORLD);
  deserialize(y, yy);
}
void mpCommunicator::join(vector<mpreal> &x, vector<mpreal> &y) {
  int M = x.size();
  int N = M/numStar;

  vector<char> xx;
  serialize(xx, x);

  int begin = my_begin*N*numDigits;
  int number = all_numbers[rank]*N*numDigits;

  vector<int> all_number = all_numbers;
  vector<int> all_begin = all_begins;

  for(int i=0; i<all_number.size(); i++) {
    all_number[i] *= numDigits*N;
    all_begin[i] *= numDigits*N;
  }

  vector<char> yy(xx.size());
  MPI_Gatherv(&xx[begin], number, MPI_CHAR, &yy.front(), &all_number.front(), &all_begin.front(), MPI_CHAR, 0, MPI_COMM_WORLD);
  deserialize(y, yy);
}

void mpCommunicator::bcast(mpreal &x, MPI_Comm comm) {
  vector<char> y(numDigits);
  convert(y, x);
  MPI_Bcast(&y.front(), numDigits, MPI_CHAR, 0, comm);  
  deconvert(x, y);
}
void mpCommunicator::bcast(vector<mpreal> &x, MPI_Comm comm) {
  int N = x.size();
  bcast(N);
  x.resize(N);
  vector<char> y(N*numDigits);
  serialize(y, x);  
  MPI_Bcast(&y.front(), N*numDigits, MPI_CHAR, 0, comm);   
  deserialize(x, y);
}
void mpCommunicator::gather(mpreal &x, vector<mpreal> &y, MPI_Comm comm) {
  y.clear();
  y.resize(size);
  vector<char> z(numDigits);
  convert(z, x);
  vector<char> col(size*numDigits);
  MPI_Gather( &z.front(), numDigits, MPI_CHAR, &col.front(), numDigits, MPI_CHAR, 0, comm);
  deserialize(y, col); 
}
void mpCommunicator::gather(vector<mpreal> &x, vector<mpreal> &y, MPI_Comm comm) {
  int N = x.size()*numDigits;
  vector<int> N_col;
  gather(N, N_col);
  vector<int> i_col(size);
  for(int i=1; i<size; i++) {
    for(int j=i; j<size; j++) {
      i_col[j] += N_col[i-1];
    }
  }
  vector<char> xx;
  serialize(xx, x);
  int N_sum = 0;
  for(int i=0; i<N_col.size(); i++) N_sum += N_col[i];
  vector<char> yy(N_sum);     
  MPI_Gatherv(&xx.front(), N, MPI_CHAR, &yy.front(), &N_col.front(), &i_col.front(), MPI_CHAR, 0, comm);
  deserialize(y, yy);
}
void mpCommunicator::join(vector<mpreal> &x, vector<mpreal> &y, MPI_Comm comm) {
  int M = x.size();
  int N = M/numStar;

  vector<char> xx;
  serialize(xx, x);

  int begin = my_begin*N*numDigits;
  int number = all_numbers[rank]*N*numDigits;

  vector<int> all_number = all_numbers;
  vector<int> all_begin = all_begins;
  for(int i=0; i<all_number.size(); i++) {
    all_number[i] *= numDigits*N;
    all_begin[i] *= numDigits*N;
  }

  vector<char> yy(xx.size());
  MPI_Gatherv(&xx[begin], number, MPI_CHAR, &yy.front(), &all_number.front(), &all_begin.front(), MPI_CHAR, 0, comm);
  deserialize(y, yy);
}
void mpCommunicator::join_by_sum(vector<mpreal> &x, vector<mpreal> &y) {
  // Make a matrix of all vectors
  vector<mpreal> data;
  gather(x, data);

  // Sum all elements
  int N = x.size();
  int p = data.size()/N;
  y.assign(N, "0");
  for(int j=0; j<p; j++) {
    for(int i=0; i<N; i++) {
      y[i] += data[j*N + i];
    }
  }
}


