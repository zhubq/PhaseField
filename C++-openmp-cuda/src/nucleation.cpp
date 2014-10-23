#include "nucleation.h"
#include <cassert>
void nucleation::Nucleate1in2(phase_field & F) {
    const int num=F.Nuc_num();
    ranlib::Uniform<float> random_number;
    Array<float,3> flag_nucleation(F.index_grain1.extent(1),F.index_grain1.extent(2),F.index_grain1.extent(3),fortranArray);
    flag_nucleation=0.;

    for(int k=F.index_grain1.lbound(3); k<=F.index_grain1.ubound(3); k++)
        for(int j=F.index_grain1.lbound(2); j<=F.index_grain1.ubound(2); j++)
            for(int i=F.index_grain1.lbound(1); i<=F.index_grain1.ubound(1); i++) {
                flag_nucleation(i,j,k)=float(count(F.index_grain1(Range::all(),i,j,k)>0))-float(F.value_conc1(2,i,j,k))*(1.-0.001*random_number.random());
            }

//         assert(any(flag_nucleation>1));
    TinyVector<int,3> Extent=flag_nucleation.extent();

    for(int n=1; n<=num; n++) {
        TinyVector<int,3> MaxIndex=maxIndex(flag_nucleation);
        float radius=1*F.I_intf();

        for(int k=1; k<=Extent(2); k++)
            for(int j=1; j<=Extent(1); j++)
                for(int i=1; i<=Extent(0); i++) {
                    float distance(TinyVector<int,3> L1, TinyVector<int,3>L2,TinyVector<int,3>Extent);
                    float Dist=distance(MaxIndex,TinyVector<int,3>(i,j,k),Extent);

                    if(Dist<radius) {
                        F.index_grain1(Range::all(),i,j,k)=0;
                        F.index_grain1(1,i,j,k)=n+100000;
                        F.value_grain1(Range::all(),i,j,k)=0.;
                        F.value_grain1(1,i,j,k)=1.;
                        F.value_phase(Range::all(),i,j,k)=0.;
                        F.value_phase(1,i,j,k)=1.;
                        //F.total_conc(i,j,k)=0.;
                        flag_nucleation(i,j,k)=0.;

                    } else if(Dist<2*radius) {
                        flag_nucleation(i,j,k)=0.;
                    }
                }
    }
}

float distance(TinyVector<int,3> L1, TinyVector<int,3>L2,TinyVector<int,3>Extent) {
    TinyVector<float,3>SqrD1=sqr(L1-L2),SqrD2=sqr(Extent-abs(L1-L2));
    float dist=0.;

    for(int i=0; i<3; i++) {
        dist+=std::min(SqrD1(i),SqrD2(i));
    }

    return sqrt(dist);
}