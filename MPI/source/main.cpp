#include <iostream>
#include "MPI_Class.h"
#include "phase_field.h"
#include "nucleation.h"
using namespace std;
int main(){
    MPI_Init(NULL,NULL);
    
    cout<<"THIS IS A PHASE FIELD SIMULATION PACKAGE (PHASPE) DEVELOPED BY BENQIANG ZHU (UNIVERSITY OF BRITISH COLUMBIA). \n"<<
          "IT IS USED FOR AT MOST THREE PHASES WHERE THE 3RD PHASE CAN BE A STOICHIOMETRIC PHASE."<<
          "FOR ANY QUESTION, PLEASE CONTACT ZHUBENQIANG@GMAIL.COM"<<endl;
    phase_field PFM;
    PFM.Read_data(("log_file1.txt"));
    int error=0;
    MPI_Barrier(MPI_COMM_WORLD);
 
    PFM.Construct(true, "dd");  
       
    PFM.Write_data(true);MPI_Barrier(MPI_COMM_WORLD);//cout<<"wsfkkds"<<endl;

    nucleation nuc;MPI_Barrier(MPI_COMM_WORLD);//cout<<"wsfdpos"<<endl;
    nuc.Nucleate1in2(PFM,1);MPI_Barrier(MPI_COMM_WORLD);//cout<<"wkkksfds"<<endl;
    PFM.Write_data(true);
    while(PFM.time<PFM.End_time())
    {
       
        PFM.Update();
        PFM.Write_data();
    }
    
    cout<<"SIMULATION IS FINISHED.\n THANK YOU FOR USING PHASPE"<<endl;

    MPI_Finalize();
    return 0;
 
}
