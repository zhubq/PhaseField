#include "omp.h"
#include <iostream>
//#include <vector>
//#include <string>
//#include <algorithm>
#include "phase_field.h"
#include "nucleation.h"
#include <time.h>
using namespace std;
int main(int argc,char**argv) {
    cout<<"THIS IS A PHASE FIELD SIMULATION PACKAGE (PHASPE) DEVELOPED BY BENQIANG ZHU (UNIVERSITY OF BRITISH COLUMBIA). \n"<<
        "IT IS USED FOR AT MOST THREE PHASES WHERE THE 3RD PHASE CAN BE A STOICHIOMETRIC PHASE."<<
        "FOR ANY QUESTION, PLEASE CONTACT ZHUBENQIANG@GMAIL.COM"<<endl;
    phase_field PFM;

    if(argc==1) {
        cout<<"PROGRAM ARGUMENT IS NOT ENOUGH."<<endl;
        exit(0);
    }

    PFM.Read_data((argv[1]));

    if(argc==2)argv[2]="dd";

    PFM.Construct(true, argv[2]);
    PFM.Write_data(true);
    nucleation nuc;
    nuc.Nucleate1in2(PFM);
    clock_t start_time=clock();
    double wstart_time=omp_get_wtime();

    while(PFM.Time()<PFM.End_time()) {
        PFM.Update(false);
        PFM.Write_data();
    }
    void GPU_Cleanup();
    cout<<"SIMULATION IS FINISHED.\n"<<endl;
    clock_t end_time=clock();
    double wend_time=omp_get_wtime();
    cout<<"Total CPU Running time(s)[physical_time*Num_threads]: "<<(end_time-start_time)/(float)CLOCKS_PER_SEC<<endl;
    cout<<"Total Physical Time (s): "<<wend_time-wstart_time<<endl;
    return 0;
}
