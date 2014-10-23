#ifndef THERMO_DATA_H
#define THERMO_DATA_H
#include <string>
#include "blitz/array.h"
#include "thermo_functor.h"


class thermo_data{
public:
    thermo_data();
    void Read_data(string);
    void Get_data(double);
    blitz::Array<double,2> para_entropy;
    blitz::Array<double,2> para_equilibrium_conc;
    blitz::Array<double,2> para_equilibrium_dc_dT;
    blitz::Array<double,2> nple1_equilibrium_conc;
    blitz::Array<double,2> nple2_equilibrium_conc;
    blitz::Array<double,2> para_partition_k;
    blitz::Array<double,2> para_partition_c;
private:
    thermo_functor functor;
};
#endif