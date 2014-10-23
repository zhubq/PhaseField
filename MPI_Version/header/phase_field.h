#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "omp.h"
 #include <assert.h>
#include "blitz/array.h"
#include "constant.h"
#include "thermo_data.h"
#include "MPI_Class.h"
using namespace blitz;
using namespace std;

static const int max_num_grains_per_phase=100000;
class phase_field{
    
public:
    phase_field();
    void Read_data(const char *);    
    void Construct(bool,char []);
    void Update();    
    void Write_data(bool forced=false);    
    double Temperature_function(double);    
    inline double End_time() const {return end_time;}
    inline TinyVector<int,3> Domain() const {return domain;}
private:
   /*this function return the Array index of the neighbor of Array(i,j,k).
 * left:-1; right:+1;lower:-1;upper:+1;forward:-1;backforward:+1.
 * for example, get Array(i-1,j,k), call neighbor(i,length=IX,left=-1);
 * get Array(i,j+2,k), call neighbor(j,length=IY,left=2);
 */
    inline int Get_neighbor(const int i,const int length,const int left) const{
//       int ghost_location=i+left;
//       if (ghost_location>length)
//         ghost_location=ghost_location-length;
//       else if (ghost_location<1)
//         ghost_location=ghost_location+length;
  
      return i+left;//ghost_location;
    }
    
    inline int Phase_of_grain(const int grain) const{
#ifdef _DEBUG
      assert(grain>0)
#endif
       return (grain/max_num_grains_per_phase>0?grain/max_num_grains_per_phase:1);
    }

    void Phase_field_solver(const int,const int,const int,const double);
    void Iso_growth(const int,const int,const int);
    void Diffusion_solver(const int,const int,const int,const double);
    double Driving_force(const int,const int,const int,const int,const double,const double) const;
    void Partitioning(const double,const Array<double,1>&,Array<double,1>&);
    void Sort_xyz(Array<double,1>& value,Array<int,1>& grain,const int N);   
    void Vti_write(ofstream &);
    void Pvti_write(ofstream & file,const int);
    void Update_parameter();
    void Tessellation(Array<int,3> & grain_index,const int num, const double reduction);
    
public:
    mpi_configuration mpi_config;
    int ix,iy,iz,i_intf;//ix,iy,iz: domain size; i_intf: number of points across interface
    int lx,ux,ly,uy,lz,uz;//subdomain bounds 
    double temperature,time,dx,eta,dimensions;//dx: grid spacing; eta: interface thickness; dt: time_step
    Array<double,2> intf_energy;
    Array<double,2> mobility,mobility0,mobility_Q;
    Array<double,1> C_diffusivity,C_diffusivity0,C_diffusivity_Q;
    //four-dimensional array for non-zero phase_field indices
    Array<int,4> index_grain1;
    Array<int,4> index_grain2;
   //four-dimensional arrayfor 
    
    Array<double,4> value_grain1;
    Array<double,4> value_grain2;
    Array<double,4> value_phase;
    Array<double,4> value_conc1;
    Array<double,4> value_conc2;
    Array<double,3> total_conc;
    enum mpi_transfer_order{INDEX,VALUE,CONC,PHASE};
    
    Array<double,1>rex_energy;
    thermo_data thermodynamic;
    double cutoff;
    string directory;
    int num_grain[max_num_phase];
    double phase_fraction[max_num_phase];
    double init_conc[max_num_phase];
    Array<double,2> temperature_profile;
    double reduction;
    int num_init_interation;
    double timestep_output;
    double end_time;
    bool is_micro_out;
    int random_seed;
    double nucleus_radius;
    char stored_energy_file[];
    int num_rex_nuclei;
    double critical_stored_energy;
    double austenite_nucleation_T;
    double num_aust_nucl_pear;
    double factor_aust_nucl_ferr;
    double pfm_dt,diff_dt;
    fstream kinetics;
    TinyVector<int,3> domain;
   
};
// void   tessellation(Array<int,3> & grain_index,int num, double reduction);
// double  distance(TinyVector<double,3>& a,TinyVector<double,3> & b,double reduction);
# endif
