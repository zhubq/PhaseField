#ifndef THERMO_FUNCTOR
#define THERMO_FUNCTOR
#include <fstream>
#include <iostream>
#include <string>
#include "blitz/array.h"
#include "constant.h"
#include "linear_functor.h"
#include "quadratic_functor.h"
#include "real.h"
using namespace std;
class thermo_functor {

public:
    thermo_functor();
    void Read_data(string);
    blitz::Array<double,2> test;
    blitz::Array<linear_functor,2> para_entropy;
    blitz::Array<quadratic_functor,2> para_equilibrium_conc;
    blitz::Array<linear_functor,2> para_equilibrium_dc_dT;
    blitz::Array<quadratic_functor,2> nple1_equilibrium_conc;
    blitz::Array<quadratic_functor,2> nple2_equilibrium_conc;
};
#endif