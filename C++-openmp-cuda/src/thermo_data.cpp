
#include "thermo_data.h"
#include <string>
#include "blitz/array.h"

using namespace blitz;
thermo_data::thermo_data():functor(),
    para_entropy(max_num_phase,max_num_phase,fortranArray),
    para_equilibrium_conc(max_num_phase,max_num_phase,fortranArray),
    para_equilibrium_dc_dT(max_num_phase,max_num_phase,fortranArray),
    nple1_equilibrium_conc(max_num_phase,max_num_phase,fortranArray),
    nple2_equilibrium_conc(max_num_phase,max_num_phase,fortranArray),
    para_partition_k(max_num_phase,max_num_phase,fortranArray),
    para_partition_c(max_num_phase,max_num_phase,fortranArray) {
    para_entropy=0.;
    para_equilibrium_conc=0.;
    para_equilibrium_dc_dT=1e-8;
}


void thermo_data::Read_data(string thermo_file) {
    functor.Read_data(thermo_file);
}
void thermo_data::Get_data(real temperature) {
    for(int j=1; j<=max_num_phase; j++)
        for(int i=1; i<=max_num_phase; i++) {
            para_entropy(i,j)=functor.para_entropy(i,j)(temperature);
            para_equilibrium_conc(i,j)=functor.para_equilibrium_conc(i,j)(temperature);
            para_equilibrium_dc_dT(i,j)=functor.para_equilibrium_dc_dT(i,j)(temperature);
            nple1_equilibrium_conc(i,j)=functor.nple1_equilibrium_conc(i,j)(temperature);
            nple2_equilibrium_conc(i,j)=functor.nple2_equilibrium_conc(i,j)(temperature);
        }

    for(int i=1; i<=max_num_phase; i++)
        for(int j=1; j<=max_num_phase; j++) {
            if(i==j) para_equilibrium_dc_dT(i,j)=1.;

            para_partition_k(i,j)=para_equilibrium_dc_dT(j,i)/para_equilibrium_dc_dT(i,j);
            para_partition_c(i,j)=para_equilibrium_conc(j,i)-para_partition_k(i,j)*para_equilibrium_conc(i,j);
        }
}