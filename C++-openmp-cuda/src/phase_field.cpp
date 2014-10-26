#include <assert.h>
#include <fstream>
#include <vector>
#include "phase_field.h"
#include <time.h>
static int GPU_Flag;
using namespace blitz;
using namespace std;
phase_field::phase_field():      intf_energy(max_num_phase,max_num_phase,fortranArray),
    mobility(max_num_phase,max_num_phase,fortranArray),
    mobility0(max_num_phase,max_num_phase,fortranArray),
    mobility_Q(max_num_phase,max_num_phase,fortranArray),
    C_diffusivity(max_num_phase,fortranArray),
    C_diffusivity0(max_num_phase,fortranArray),
    C_diffusivity_Q(max_num_phase,fortranArray),
    index_grain1(fortranArray),
    index_grain2(fortranArray),
    value_grain1(fortranArray),
    value_grain2(fortranArray),
    value_phase(fortranArray),
    value_conc1(fortranArray),
    value_conc2(fortranArray),
    total_conc(fortranArray),
    temperature_profile(fortranArray),
    rex_energy(fortranArray) {
    time=0.;
    pfm_dt=0.;
    temperature=0;
}

void phase_field::Read_data(const char filename[]="read_from_keyboard") {
    GPU_Flag=0;
    string thermo_file;

    if (filename=="read_from_keyboard") {
        cin>>ix>>iy>>iz;
        cin>>i_intf;
        cin>>dx;
        cin>>thermo_file;

    } else {
        cout<<"file name is "<<filename<<endl;
        fstream input_file;
        input_file.open(filename);
        input_file.ignore(1000,'\n');
        getline(input_file,directory,'\n');
        string command="mkdir ";
        const char* a=(command+directory).c_str();
        cout<<"directory is "<<directory<<'.'<<endl;
        system(a);
        kinetics.open((directory+"/fraction_time.txt").c_str(),ios::out);
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>ix>>iy>>iz;
        cout<<"domain size is: "<<ix<<" "<<iy<<" "<<iz<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>dx;
        cout<<"cell size (um): "<<dx<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>i_intf;
        cout<<"number of grid points across the interface"<<i_intf<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>cutoff;
        cout<<"cutoff:"<<cutoff<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"grain numbers in the initial microstructure:"<<endl;

        for(int i=0; i<max_num_phase; i++) {
            input_file>>num_grain[i];
            cout<<num_grain[i]<<" ";
        }

        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"\nPhase fractions in the initial microstructure:"<<endl;

        for(int i=0; i<max_num_phase; i++) {
            input_file>>phase_fraction[i];
            cout<<phase_fraction[i]<<" ";
        }

        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"\nCarbon concentrations in the initial microstructure:"<<endl;

        for(int i=0; i<max_num_phase; i++) {
            input_file>>init_conc[i];
            cout<<init_conc[i]<<" ";
        }

        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        int num_segment_T;
        input_file>>num_segment_T;
        temperature_profile.resize(2,num_segment_T);

        for(int i=1; i<=num_segment_T; i++) {
            input_file>>temperature_profile(1,i)>>temperature_profile(2,i);
            cout<<temperature_profile(1,i)<<" "<<temperature_profile(2,i)<<endl;
        }

        time=temperature_profile(1,1);
        temperature=temperature_profile(2,1);
        end_time=temperature_profile(1,num_segment_T);
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>reduction;
        cout<<"\ncold rolling reduction:"<<reduction<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>num_init_interation;
        cout<<"number of initial interation to generate diffuse interface is "
            <<num_init_interation<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>timestep_output;
        cout<<"time interval to output simulation results: "<<timestep_output<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>Is_micro_out;

        if(Is_micro_out)cout<<"microstructure files(vti) are output."<<endl;

        input_file.ignore(1000,'\n');
//         input_file.ignore(1000,'\n');
//         input_file>>random_seed;
//         input_file.ignore(1000,'\n');
//         input_file.ignore(1000,'\n');
//         input_file>>nucleus_radius;
//         input_file.ignore(1000,'\n');
//         input_file.ignore(1000,'\n');
//         input_file.getline(stored_energy_file,'\n');
//         cout<<stored_energy_file<<'.'<<endl;
//         input_file.ignore(1000,'\n');
//         input_file>>num_rex_nuclei;
//         cout<<"number of ferrite recrystallization nuclei: "<<num_rex_nuclei<<endl;
//         input_file.ignore(1000,'\n');
//         input_file.ignore(1000,'\n');
//         input_file>>critical_stored_energy;
//         cout<<"critical stored energy for recrystallization nucleation is: "
//             <<critical_stored_energy<<" MPa"<<endl;
//         input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
//         input_file>>austenite_nucleation_T>>num_aust_nucl_pear>>factor_aust_nucl_ferr;
        input_file>>num_aust_nucl_ferr;

        for(int i=1; i<=max_num_phase; i++)
            for(int j=i; j<=max_num_phase; j++) {
                input_file.ignore(1000,'\n');
                input_file.ignore(1000,'\n');
                input_file>>mobility0(i,j)>>mobility_Q(i,j);
                mobility0(j,i)=mobility0(i,j);
                mobility_Q(j,i)=mobility_Q(i,j);
                cout<<"interface energy between phase "<<i<<" and "<<j<<" is "
                    <<mobility0(i,j)<<"um4/J/s and "<<mobility_Q(i,j)<<"J/mol"<<endl;
            }

        for(int i=1; i<=max_num_phase; i++)
            for(int j=i; j<=max_num_phase; j++) {
                input_file.ignore(1000,'\n');
                input_file.ignore(1000,'\n');
                input_file>>intf_energy(i,j);
                intf_energy(j,i)=intf_energy(i,j);
                cout<<"interface energy between phase "<<i<<" and "<<j<<" is "
                    <<intf_energy(i,j)<<"J/um2"<<endl;
            }

        for(int i=1; i<=max_num_phase; i++) {
            input_file.ignore(1000,'\n');
            input_file.ignore(1000,'\n');
            input_file>>C_diffusivity0(i)>>C_diffusivity_Q(i);
            cout<<"carbon diffusivity in phase "<<i<<" : "<<C_diffusivity0(i)<<"um2/s "
                <<C_diffusivity_Q(i)<<"J/mol"<<endl;
        }

        for(int i=1; i<=2; i++)input_file.ignore(1000,'\n');

        //char thermo_file[];
        getline(input_file,thermo_file,'\n');
        cout<<"thermo-data file is: "<<thermo_file.c_str()<<'.'<<endl;
        input_file.ignore(1000,'\n');
        //input_file.ignore(1000,'\n');
        input_file>>GPU_Flag;

        if(GPU_Flag)cout<<"GPU will be used if available"<<endl;

        //exit(0);
        ;
    }

    thermodynamic.Read_data(thermo_file);
    index_grain1.resize(max_num_coexist,ix,iy,iz);
    index_grain1=0.;
    index_grain2.resize(max_num_coexist,ix,iy,iz);
    index_grain2=0;
    value_grain1.resize(max_num_coexist,ix,iy,iz);
    value_grain1=0.;
    value_grain2.resize(max_num_coexist,ix,iy,iz);
    value_grain2=0.;
    value_phase.resize(max_num_phase,ix,iy,iz);
    value_phase=0.;
    value_conc1.resize(max_num_phase,ix,iy,iz);

    for(int i=1; i<=max_num_phase; i++) {
        value_conc1(i,Range::all(),Range::all(),Range::all())=init_conc[i-1];
    }

    value_conc2.resize(max_num_phase,ix,iy,iz);
    value_conc2=0.;
    total_conc.resize(ix,iy,iz);
    total_conc=0.;
    rex_energy.resize(num_grain[0]);
    rex_energy=0.;
    dimensions=3;

    if(iz==1)
        dimensions-=1;

    if(iy==1)
        dimensions-=1;

    cout<<"the dimensionality of the simulation domain is "<<dimensions<<endl;
    eta=dx*i_intf;
    GPU_Flag=GPU_Flag&&GPU_Init(ix,iy,iz);
    return;
}
real phase_field::Temperature_function(real Time) {
    real Temp;

    for(int i=1; i<temperature_profile.extent(1); i++)
        if(Time>=temperature_profile(1,i)&&Time<temperature_profile(1,i+1)) {
            Temp=(temperature_profile(2,i+1)-temperature_profile(2,i))/(temperature_profile(1,i+1)-temperature_profile(1,i))
                 *(Time-temperature_profile(1,i))+temperature_profile(2,i);
        }

    return Temp;
}
void phase_field::Update(bool flag_diff) {
    int num_iter_diff_per_pfm;
    real new_diff_dt;
    phase_field::Update_parameter();
    num_iter_diff_per_pfm=max(int(ceil(pfm_dt/diff_dt)),1);
    new_diff_dt=pfm_dt/num_iter_diff_per_pfm;
//     clock_t s1=clock();
//     double ss1=omp_get_wtime();
//USE FINITE-DIFFERENCE TO SOLVE PHASE FIELD EQUATIONS ON EACH GRID POINT.
    #pragma omp parallel for schedule(guided,1) collapse(2)

    for(int k=1; k<=iz; k++)
        for(int j=1; j<=iy; j++)
            for(int i=1; i<=ix; i++)
                Phase_field_solver(i,j,k,pfm_dt);

    cycleArrays(index_grain1,index_grain2);//SWITCH POINTERS
    cycleArrays(value_grain1,value_grain2);

    if(flag_diff)return;

//     clock_t s2=clock();
//     double ss2=omp_get_wtime();
//     cout<<"pfm "<<ss2-ss1<<endl;
    if(!GPU_Flag) {
        /*************************************************
         * *****CPU DIFFUSION SOLVER**********************
         */
	
//#pragma omp parallel
        for(int n=1; n<=num_iter_diff_per_pfm; n++) {
            #pragma omp parallel for schedule(static) collapse(2)

            for(int k=1; k<=iz; k++)
                for(int j=1; j<=iy; j++)
                    for(int i=1; i<=ix; i++)
                        Diffusion_solver(i,j,k,new_diff_dt);

//#pragma omp single
            cycleArrays(value_conc1,value_conc2);//cout<<total_conc(1,1,1)<<" "<<value_conc1(2,1,1,1)<<" "<<value_phase(2,1,1,1)<<"/";
        }

    } else {
        /*************************************************
         *********GPU DIFFUSION SOLVER********************
         */
        GPU_Diffusion_solver(C_diffusivity.data(), thermodynamic.para_partition_k.data(),thermodynamic.para_partition_c.data(),
                             value_phase.data(),total_conc.data(),ix,iy, iz,pfm_dt,dx);
        #pragma omp for collapse(3)

        for(int k=1; k<=iz; k++)
            for(int j=1; j<=iy; j++)
                for(int i=1; i<=ix; i++) {
                    Partitioning(total_conc(i,j,k),&value_phase(1,i,j,k),&value_conc1(1,i,j,k));
                }
    }

    //return;
//     clock_t s3=clock();
//     double ss3=omp_get_wtime();
//     cout<<"diff "<<ss3-ss2<<num_iter_diff_per_pfm<<endl;
    return;
}

void phase_field::Update_parameter() {
    temperature=Temperature_function(time);
    mobility=mobility0*exp(-mobility_Q/(R*temperature));
    C_diffusivity=C_diffusivity0*exp(-C_diffusivity_Q/(R*temperature));
    pfm_dt=dx*dx/dimensions*0.2/(blitz::max(mobility*intf_energy));
    diff_dt=dx*dx/blitz::max(C_diffusivity)/dimensions*0.2;
    thermodynamic.Get_data(temperature);
    time+=pfm_dt;
}

