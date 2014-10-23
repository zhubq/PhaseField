#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "omp.h"
// #include <mpi.h>
#include "blitz/array.h"
#include "constant.h"
#include "thermo_data.h"
#include "gpu_solver.h"
#include "real.h"
using namespace blitz;
using namespace std;

static const int max_num_grains_per_phase=100000;
class phase_field {
public:
    phase_field();//constructor;
    void Read_data(const char *);//read in data from input files and initialize pointers;    
    //get the phase index for a grain index
    inline int Phase_of_grain(const int grain) const {
        if(grain==0)return 0;

        return (grain/max_num_grains_per_phase>0?grain/max_num_grains_per_phase:1);
    }
    inline const real Time() const {
        return time;
    }

    inline const real Temperature() const {
        return temperature;
    }

    inline const real End_time() const {
        return end_time;
    }
    inline const int I_intf() const {
        return i_intf;
    }
    inline const int Nuc_num() const {
        return num_aust_nucl_ferr;
    }

    void Construct(bool,char []);//construct initial microstructure;
    void Update(bool);//update phase field and solute concentration for each time step;
    void Write_data(bool forced=false);//write results to files for postprocessing;

private:
      /*this function return the Array index of the neighbor of Array(i,j,k).
    * left:-1; right:+1;lower:-1;upper:+1;forward:-1;backforward:+1.
    * for example, get Array(i-1,j,k), call neighbor(i,length=IX,left=-1);
    * get Array(i,j+2,k), call neighbor(j,length=IY,left=2);
    */
    inline int Get_neighbor(const int i,const int length,const int left) const {
        int ghost_location=i+left;

        if (ghost_location>length)
            ghost_location=ghost_location-length;

        else if (ghost_location<1)
            ghost_location=ghost_location+length;

        return ghost_location;
    }
    void Tessellation(Array<int,3> & grain_index,int num, real reduction);//construct a polycrystal structure using Voronoi tessellation;
    void Iso_growth(const int,const int,const int);//Generate diffuse interfaces assuming isotropic growth of initial polycrystal grains;
    void Diffusion_solver(const int i,const int j,const int k,const real time_step);//solve diffusion equations in all phases on each grid point(i,j,k);
    void Phase_field_solver(const int i,const int j,const int k,const real time_step);//solve phase-field equations on each grid point (i,j,k);
    void Partitioning(const real,const real*,real*); //partition total concentrations to each phase concentration;  
    //void Partitioning(const real,const Array<real,1>&,Array<real,1>&);
    real Driving_force(const int index_ii,const int phase_ii,const int index_jj,const int phase_jj,const real conc_ii,const real conc_jj) const;//get the driving pressure of two interacting phase fields;    
    void Sort_xyz(Array<real,1>& value,Array<int,1>& grain,const int N);//sort a vector of phase-field values and relevant indices;
    void Update_parameter();//within each time step, update all relevant parameters, e.g.time, temperature,diffusivity,mobility and thermodynamic data;
    real Temperature_function(real);//return the temperature given a time, based on the inputted temperature path;
    void Vti_write(ofstream &,Array<int,4> &,Array<real,3> &,real);// ouput VTI-formatted files, both the polycrystal structure and the concentration field;
public:    
    Array<int,4> index_grain1; //4D matrix, stores non-zero phase-field indices on each 3D grid points;
    Array<real,4> value_grain1; //4D matrix, stores values of non-zero indices on each 3D grid points;
    Array<real,4> value_phase; // 4D matrix, stores fractions of phases on each 3D grid points;
    Array<real,4> value_conc1; //4D matrix, stores solute concentration of phases on each 3D grid points;
    Array<real,3> total_conc; //3D matrix, stores total concentration on each 3D grid points, equals dot_product(value_phase,value_conc1);
private: 
    Array<int,4> index_grain2;   //temporary matrix, storing updated non-zero indices in each iteration;
    Array<real,4> value_grain2;   //temporary matrix; 
    Array<real,4> value_conc2;   //temporary matrix;
    Array<real,1>rex_energy;  //stored energy;

    int ix,iy,iz,i_intf;//ix,iy,iz: domain size; i_intf: number of points across interface
    real temperature,time,dx,eta,dimensions;//dx: grid spacing; eta: interface thickness; dt: time_step
    Array<real,2> intf_energy;//interface energy 2D-matrix[max_num_phase][max_num_phase]
    Array<real,2> temperature_profile;//temperature path 2D-matrix[num_of_T_points][2]: first-column--time, second-column--temperature;
    thermo_data thermodynamic;//thermodynamic data (a "thermo_data" object), include equilibrium concentrations and slopes of solidus and liquidus lines;
    Array<real,2> mobility,mobility0,mobility_Q;//interface mobilities, pre-factors and activation energies,2D-matrix [max_num_phase][max_num_phase];
    Array<real,1> C_diffusivity,C_diffusivity0,C_diffusivity_Q;//solute diffusivities, pre-factors and activation energies,1D array [max_num_phase];
    real reduction;//cold-rolling thickness reduction, used as an input parameter for tessellation to generate an elongated-grain morphology;
    int num_init_interation;//number of iterations for isogrowth to generate diffuse interfaces for the inital microstructure;
    real timestep_output;//time interval to write files;
    real end_time;//simulation time;
    bool Is_micro_out;//flag for whether write VTI files (if =1) or only statistical data (if =0);
    //int random_seed;
    //real nucleus_radius;
    //char stored_energy_file[];
    //int num_rex_nuclei;
    //real critical_stored_energy;
    //real austenite_nucleation_T;
    real num_aust_nucl_ferr; //number of nuclei;
    //real factor_aust_nucl_ferr;
    real pfm_dt,diff_dt;//time steps for phase-field and diffusion equations;
    fstream kinetics;//file stream object for output of statistical data
    real cutoff;//threshod value for interface edges: if all phase fields on a grid point is smaller than it, the grid point is meant out of the interface region;
    string directory;// folder name to store output files;
    int num_grain[max_num_phase];//number of grains in the initial microstructure, used to construct the initial microstructure;
    real phase_fraction[max_num_phase];// volume fraction of each phase in the initial microstructure, used for construction;
    real init_conc[max_num_phase];// initial concentrations of each phase,assuming uniform in each phase;

};

# endif
