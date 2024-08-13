#include "mpi.h"

#include <iostream>
using namespace std;
#include <vector>
#include <string>

#include "mpreal.h"
using namespace mpfr;

#ifndef __mpCommunicator_h
#define __mpCommunicator_h

class mpCommunicator {

  int rank, size;
  string name;

  int numStar;
  int my_begin, my_end;
  vector<int> all_begins, all_numbers;

  int numDigits;

  int group, sendbuf_rank, recvbuf_sum;
  MPI_Group orig_group, new_group1, new_group2;  
  MPI_Comm new_comm, new_comm1, new_comm2;

  public:

  // constructors
  mpCommunicator();

  // mpi
  void start_mpi(int argc, char* argv[]);
  void stop_mpi();

  void make_two_groups();
  void free_groups();

  // work
  void divide_work(int numElements);

  // mpreal
  void set_mpreal(int numDigits);
  void format(string &x);
  void convert(vector<char> &x, mpreal &y);
  void deconvert(mpreal &x, vector<char> &y);  
  void serialize(vector<char> &y, vector<mpreal> &x);
  void deserialize(vector<mpreal> &x, vector<char> &y);

  // getters
  string get_name();

  int get_rank();
  int get_size();
  int get_my_begin();
  int get_my_end();

  int get_group();
  MPI_Comm get_comm_world();
  MPI_Comm get_comm_group1();
  MPI_Comm get_comm_group2();

  // int
  void bcast(int &x);
  void bcast(vector<int> &x);
  void gather(int &x, vector<int> &y);
  void gather(vector<int> &x, vector<int> &y);
  void join(vector<int> &x, vector<int> &y);
  void send_id();
  void recv_id(int &src_id);
  void send_ic(int &dest_id, int &ic_index);
  void recv_ic(int &ic_index);

  void bcast(int &x, MPI_Comm comm);
  void bcast(vector<int> &x, MPI_Comm comm);
  void gather(int &x, vector<int> &y, MPI_Comm comm);
  void gather(vector<int> &x, vector<int> &y, MPI_Comm comm);
  void join(vector<int> &x, vector<int> &y, MPI_Comm comm);
  void send_id(MPI_Comm comm);
  void recv_id(int &src_id, MPI_Comm comm);
  void send_ic(int &dest_id, int &ic_index, MPI_Comm comm);
  void recv_ic(int &ic_index, MPI_Comm comm);

  // float
  void bcast(float &x);
  void bcast(vector<float> &x);
  void gather(float &x, vector<float> &y);
  void gather(vector<float> &x, vector<float> &y);
  void join(vector<float> &x, vector<float> &y);

  void bcast(float &x, MPI_Comm comm);
  void bcast(vector<float> &x, MPI_Comm comm);
  void gather(float &x, vector<float> &y, MPI_Comm comm);
  void gather(vector<float> &x, vector<float> &y, MPI_Comm comm);
  void join(vector<float> &x, vector<float> &y, MPI_Comm comm);

  // double
  void bcast(double &x);
  void bcast(vector<double> &x);
  void gather(double &x, vector<double> &y);
  void gather(vector<double> &x, vector<double> &y);
  void join(vector<double> &x, vector<double> &y);

  void bcast(double &x, MPI_Comm comm);
  void bcast(vector<double> &x, MPI_Comm comm);
  void gather(double &x, vector<double> &y, MPI_Comm comm);
  void gather(vector<double> &x, vector<double> &y, MPI_Comm comm);
  void join(vector<double> &x, vector<double> &y, MPI_Comm comm);

  // mpreal
  void bcast(mpreal &x);
  void bcast(vector<mpreal> &x);
  void gather(mpreal &x, vector<mpreal> &y);
  void gather(vector<mpreal> &x, vector<mpreal> &y);
  void join(vector<mpreal> &x, vector<mpreal> &y);
  void join_by_sum(vector<mpreal> &x, vector<mpreal> &y);

  void bcast(mpreal &x, MPI_Comm comm);
  void bcast(vector<mpreal> &x, MPI_Comm comm);
  void gather(mpreal &x, vector<mpreal> &y, MPI_Comm comm);
  void gather(vector<mpreal> &x, vector<mpreal> &y, MPI_Comm comm);
  void join(vector<mpreal> &x, vector<mpreal> &y, MPI_Comm comm);
};

#endif


