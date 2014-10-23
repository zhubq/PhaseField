#ifndef __MPI_CLASS
#define __MPI_CLASS
//#define _DEBUG
#include <iostream>
#include <mpi.h>
#include "blitz/array.h"
using namespace blitz;
class mpi_configuration
{
public:
  mpi_configuration();
  
  // configure 3D cartensian topology for mpi
  // decompose the domain into subdomains
  void Configure_cart(const TinyVector<int,3> &);
  // create and commit data types that are used to transfer ghost info between each pair of neighbors 
  void Commit_cart_datatype(const int,const int);
  // start non-blocking communication for ghost cells of indices of grains [index_grain1]
  void Start_nonblock_index( Array<int,4> & );
  // start non-blocking communication for ghost cells of grain fractions [value_grain1]
  void Start_nonblock_value( Array<double,4> & );
  // start non-blocking communication for ghost cells of solute concentration [total_conc]
  void Start_nonblock_conc( Array<double,3> & );
  // start non-blocking communication for ghost cells of phase fractions [value_phase]
  void Start_nonblock_phase( Array<double,4> & );
  // wait for non-blocking communications above finished
  // i represent each communication
  // 0: index  1:value  2: conc  3: phase
  void Wait_nonblock_request(int i);
  // free all nonblocking requests
  void Free_nonblock_request();
  
  inline int Cart_rank() const {return myrank;}
  inline int Cart_size() const {return comm_size;}
  inline MPI_Comm Cart_comm() const {return cart_comm;}
  inline TinyVector<int,6> Subdomain_size(const int rank) const {return subdomain_info(rank);}
//private:  
  int myrank;  //mpi rank
  int comm_size;  //mpi comm size
  MPI_Comm cart_comm;  // cartensian communicator handler
  int cart_dims[3];//Number of processes per dimension.
  int cart_periods[3]; // flag of periods: 1--periodic 0--non-periodic
  int cart_coord[3];  // coordinates of the process in cartensian topology
  int cart_neighbors[6];  // ranks of 6 neighbors in cartensian topology
  Array<TinyVector<int,6>,1> subdomain_info; // arrays storing subdomain sizes of all processes Tinyvector--[lx ux ly uy lz uz] lower and upper bounds of each direction
  // Variables for nonblocking communication
  // 3D subdomain: thus 6 faces, transfer 4 variables: index_grain, value_grain and total_conc,value_phase
  void* send_buffer_ptr[4][6]; //send buffer address
  void* recv_buffer_ptr[4][6]; //recv buffer address
  //int comm_count[4][6];
  MPI_Datatype data_type[4][6]; //datatype of ghost cells
  int nonblock_send_tag[4][6]; //
  int nonblock_recv_tag[4][6];
  MPI_Request nonblock_send_request[4][6];
  MPI_Request nonblock_recv_request[4][6];
  MPI_Status nonblock_send_status[4][6];
  MPI_Status nonblock_recv_status[4][6];
  
};

#endif