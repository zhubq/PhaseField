
#include "thermo_functor.h"
#include <string>
#include "blitz/array.h"
using namespace blitz;
thermo_functor::thermo_functor():para_entropy(max_num_phase,max_num_phase,blitz::fortranArray),
    para_equilibrium_conc(max_num_phase,max_num_phase,blitz::fortranArray),
    para_equilibrium_dc_dT(max_num_phase,max_num_phase,blitz::fortranArray),
    nple1_equilibrium_conc(max_num_phase,max_num_phase,blitz::fortranArray),
    nple2_equilibrium_conc(max_num_phase,max_num_phase,blitz::fortranArray) {}
void thermo_functor::Read_data(string thermo_file) {
    std::ifstream file;
    file.open(thermo_file.c_str());
    file.ignore(1000,'\n');
    file.ignore(1000,'\n');
    file.ignore(1000,'\n');
    int i,j;
    real a,b,c;

    for(int k=1; k<=(max_num_phase-1)*max_num_phase/2; k++) {
        file>>i>>j>>a>>b;
        //std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<std::endl;
        para_entropy(i,j).Reset(a,b);
        para_entropy(j,i).Reset(-a,-b);
    }

    for(int i=1; i<=4; i++)       file.ignore(1000,'\n');

    for(int k=1; k<=(max_num_phase-1)*max_num_phase; k++) {
        file>>i>>j>>a>>b>>c;
        //std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
        nple1_equilibrium_conc(i,j).Reset(a,b,c);
    }

    for(int i=1; i<=2; i++)file.ignore(1000,'\n');

    for(int k=1; k<=(max_num_phase-1)*max_num_phase; k++) {
        file>>i>>j>>a>>b>>c;
        //std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
        nple2_equilibrium_conc(i,j).Reset(a,b,c);
    }

    for(int i=1; i<=2; i++)file.ignore(1000,'\n');

    for(int k=1; k<=(max_num_phase-1)*max_num_phase; k++) {
        file>>i>>j>>a>>b>>c;
        //std::cout<<i<<" "<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
        para_equilibrium_conc(i,j).Reset(a,b,c);
        para_equilibrium_dc_dT(i,j).Reset(2*a,b);
    }

    return;
}