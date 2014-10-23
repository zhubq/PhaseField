#include <iostream>
#include "blitz/array.h"
#include "gpu_solver.h"
#include "real.h"
using namespace std;
using namespace blitz;

int main() {
    cout<<"hello GPU"<<endl;
    GPU_Init(320,320,1);
    real diff[3]= {1,1,1};
    Array<real,2> partk(3,3,fortranArray);
    partk=1.f,2.f,1.f,0.5f,1.f,1.f,1.f,1.f,1.f;
    real partc[3][3]= {0.f};
    //cout<<sizeof(partc)<<endl;
    Array<real,4>value_phase(3,320,320,1,fortranArray);
    Array<real,3>total_conc(320,320,1,fortranArray);
    value_phase=0.f;
    value_phase(1,Range::all(),Range::all(),Range::all())=1.0f;
    total_conc=0.1f;
    total_conc(Range::all(),Range::all(),Range::all())=0.2f;
    cout<<total_conc(12,12,Range::all())<<endl;
    //total_conc=0.f;
    real n;
    cin>>n;
    Diffusion_solver(diff,partk.data(),partc[0],value_phase.data(),total_conc.data(),320,320,1,n,1.);
    //Array<real,4>val_conc(3,320,320,1,fortranArray);
    GPU_Cleanup();
    return 0;
}