#include <assert.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "phase_field.h"

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
                                  rex_energy(fortranArray)//,
                                 // mpi_config(Config)
 {
     time=0.;
     pfm_dt=0.;
     temperature=0;
     
}

void phase_field::Read_data(const char filename[]="read_from_keyboard")
{
    string thermo_file;   
    if (filename=="read_from_keyboard") {
        cin>>ix>>iy>>iz;
        cin>>i_intf;
        cin>>dx;
        cin>>thermo_file;
    }
    else{
        cout<<"file name is "<<filename<<endl;
        fstream input_file;
        input_file.open(filename);
        input_file.ignore(1000,'\n');
        getline(input_file,directory,'\n');
        string command="mkdir ";
        const char* a=(command+directory).c_str();
        cout<<"directory is "<<directory<<'.'<<endl;
        system(a);
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        if(myrank==0)kinetics.open((directory+"/fraction_time.txt").c_str(),ios::out);
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>ix>>iy>>iz;
        domain[0]=ix;
        domain[1]=iy;
        domain[2]=iz;
        cout<<"domain size is: "<<ix<<" "<<iy<<" "<<iz<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>dx;
        cout<<"cell size (um): "<<dx<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>i_intf;
        cout<<"number of grid points across the interface"<<i_intf<<endl;
	eta=dx*i_intf;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>cutoff;cout<<"cutoff:"<<cutoff<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"grain numbers in the initial microstructure:"<<endl;
        for(int i=0;i<max_num_phase;i++) {
            input_file>>num_grain[i]; 
            cout<<num_grain[i]<<" ";
            
        }
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"\nPhase fractions in the initial microstructure:"<<endl;
        for(int i=0;i<max_num_phase;i++) {
            input_file>>phase_fraction[i];
            cout<<phase_fraction[i]<<" ";
      
        }
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        cout<<"\nCarbon concentrations in the initial microstructure:"<<endl;
        for(int i=0;i<max_num_phase;i++) {
            input_file>>init_conc[i];
            cout<<init_conc[i]<<" ";
            
        }
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        int num_segment_T; 
        input_file>>num_segment_T;
        temperature_profile.resize(2,num_segment_T);
        for(int i=1; i<=num_segment_T;i++){
            
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
        input_file>>is_micro_out;
        if(is_micro_out)cout<<"microstructure files(vti) are output."<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>random_seed;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>nucleus_radius;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file.getline(stored_energy_file,'\n');
        cout<<stored_energy_file<<'.'<<endl;
        input_file.ignore(1000,'\n');//input_file.ignore(1000,'\n');
        input_file>>num_rex_nuclei;
        cout<<"number of ferrite recrystallization nuclei: "<<num_rex_nuclei<<endl; 
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>critical_stored_energy;
        cout<<"critical stored energy for recrystallization nucleation is: "
            <<critical_stored_energy<<" MPa"<<endl;
        input_file.ignore(1000,'\n');
        input_file.ignore(1000,'\n');
        input_file>>austenite_nucleation_T>>num_aust_nucl_pear>>factor_aust_nucl_ferr;
        for(int i=1;i<=max_num_phase;i++)
          for(int j=i;j<=max_num_phase;j++){
            input_file.ignore(1000,'\n');input_file.ignore(1000,'\n');
            input_file>>mobility0(i,j)>>mobility_Q(i,j);
            mobility0(j,i)=mobility0(i,j);
            mobility_Q(j,i)=mobility_Q(i,j);
            cout<<"interface energy between phase "<<i<<" and "<<j<<" is "
                <<mobility0(i,j)<<"um4/J/s and "<<mobility_Q(i,j)<<"J/mol"<<endl;
        }
        for(int i=1;i<=max_num_phase;i++)
          for(int j=i;j<=max_num_phase;j++){
            input_file.ignore(1000,'\n');input_file.ignore(1000,'\n');
            input_file>>intf_energy(i,j);
            intf_energy(j,i)=intf_energy(i,j);
            cout<<"interface energy between phase "<<i<<" and "<<j<<" is "
                <<intf_energy(i,j)<<"J/um2"<<endl;
        }
        for(int i=1;i<=max_num_phase;i++){
          input_file.ignore(1000,'\n');
          input_file.ignore(1000,'\n');
          input_file>>C_diffusivity0(i)>>C_diffusivity_Q(i);
          cout<<"carbon diffusivity in phase "<<i<<" : "<<C_diffusivity0(i)<<"um2/s "
              <<C_diffusivity_Q(i)<<"J/mol"<<endl;
        }
        for(int i=1;i<=2;i++)input_file.ignore(1000,'\n');
        //char thermo_file[];
        
        getline(input_file,thermo_file,'\n');
        
        cout<<"thermo-data file is: "<<thermo_file.c_str()<<'.'<<endl;
        //exit(0);
      ;
    }
    
    thermodynamic.Read_data(thermo_file);
    
    dimensions=3;
    if(iz==1)
      dimensions-=1;
    if(iy==1)
      dimensions-=1;
    cout<<"the dimensionality of the simulation domain is "<<dimensions<<endl;
    cout<<"domain size "<<domain<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_config.Configure_cart(domain);//cout<<"1"<<endl;
    //MPI_Barrier(mpi_config.Cart_comm());cout<<mpi_config.myrank<<mpi_config.cart_neighbors[0]<<mpi_config.cart_neighbors[1]<<mpi_config.cart_neighbors[2]
   // <<mpi_config.cart_neighbors[3]<<mpi_config.cart_neighbors[4]<<mpi_config.cart_neighbors[5]<<endl;
    mpi_config.Commit_cart_datatype(max_num_coexist,max_num_phase);
    MPI_Barrier(mpi_config.Cart_comm());//cout<<3<<endl;
    TinyVector<int,6> mysubdomain=mpi_config.Subdomain_size(mpi_config.Cart_rank());
    lx=mysubdomain(0);
    ux=mysubdomain(1);
    ly=mysubdomain(2);
    uy=mysubdomain(3);
    lz=mysubdomain(4);
    uz=mysubdomain(5);
    
    TinyVector<int,4> baseIndex(1,lx,ly,lz);
    index_grain1.reindexSelf(baseIndex);
    index_grain2.reindexSelf(baseIndex);
    value_grain1.reindexSelf(baseIndex);
    value_grain2.reindexSelf(baseIndex);
    value_conc1.reindexSelf(baseIndex);
    value_conc2.reindexSelf(baseIndex);
    total_conc.reindexSelf(TinyVector<int,3>(lx,ly,lz));
    value_phase.reindexSelf(baseIndex);
    index_grain1.resize(max_num_coexist,ux-lx+1,uy-ly+1,uz-lz+1);index_grain1=0.;
    index_grain2.resize(max_num_coexist,ux-lx+1,uy-ly+1,uz-lz+1);index_grain2=0;
    value_grain1.resize(max_num_coexist,ux-lx+1,uy-ly+1,uz-lz+1);value_grain1=0.;
    value_grain2.resize(max_num_coexist,ux-lx+1,uy-ly+1,uz-lz+1);value_grain2=0.;
    value_phase.resize(max_num_phase,ux-lx+1,uy-ly+1,uz-lz+1);value_phase=0.;
    value_conc1.resize(max_num_phase,ux-lx+1,uy-ly+1,uz-lz+1);
    for(int i=1;i<=max_num_phase;i++){
        value_conc1(i,Range::all(),Range::all(),Range::all())=init_conc[i-1];
    }
    
    value_conc2.resize(max_num_phase,ux-lx+1,uy-ly+1,uz-lz+1);value_conc2=0.;
    total_conc.resize(ux-lx+1,uy-ly+1,uz-lz+1);total_conc=0.;
    //rex_energy.resize(num_grain[0]);rex_energy=0.;

    assert(index_grain1.isStorageContiguous());
    assert(index_grain2.isStorageContiguous());
    assert(value_grain1.isStorageContiguous());
    assert(value_grain2.isStorageContiguous());
    assert(value_phase.isStorageContiguous());
    assert(value_conc1.isStorageContiguous());
    assert(value_conc2.isStorageContiguous());
    
    if(mpi_config.Cart_rank()==0)cout<<"phase_field.Read_data() is successfully called."<<endl;
    
   

    return;
}
double phase_field::Temperature_function(double Time)
{
  double Temp;
  for(int i=1;i<temperature_profile.extent(1);i++)
    if(Time>=temperature_profile(1,i)&&Time<temperature_profile(1,i+1)){
      Temp=(temperature_profile(2,i+1)-temperature_profile(2,i))/(temperature_profile(1,i+1)-temperature_profile(1,i))
                 *(Time-temperature_profile(1,i))+temperature_profile(2,i);
		 
    }
  return Temp;  
}
void phase_field::Update()
{
  
  int n;
  int num_iter_diff_per_pfm;
  double new_diff_dt;
  phase_field::Update_parameter();
  num_iter_diff_per_pfm=max(int(ceil(pfm_dt/diff_dt)),1);
  new_diff_dt=pfm_dt/num_iter_diff_per_pfm;
  //if(mpi_config.Cart_rank()==0)cout<<time<<" "<<num_iter_diff_per_pfm<<endl;
  mpi_config.Start_nonblock_index(index_grain1);
  mpi_config.Start_nonblock_value(value_grain1);
      
  #pragma omp parallel for collapse(3)
      // calculate on the core points
  for(int k=lz+2;k<=uz-2;k++)
     for(int j=ly+2;j<=uy-2;j++)
        for(int i=lx+2;i<=ux-2;i++)
          {
              
              Phase_field_solver(i,j,k,pfm_dt);
          }
  mpi_config.Wait_nonblock_request(INDEX);
  mpi_config.Wait_nonblock_request(VALUE);
      
 
      // calculate surface points

      //XY PLANES
  for(int j=ly+1;j<=uy-1;j++)
          for(int i=lx+1;i<=ux-1;i++)
          {   //cout<<"clc ";
            Phase_field_solver(i,j,lz+1,pfm_dt);
            if(uz-lz>2)Phase_field_solver(i,j,uz-1,pfm_dt);
          } 
      //YZ PLANES
  for(int k=lz+2;k<=uz-2;k++)
          for(int j=ly+1;j<=uy-1;j++)
          {
            Phase_field_solver(lx+1,j,k,pfm_dt);
            if(ux-lx>2) Phase_field_solver(ux-1,j,k,pfm_dt);
          }
      //ZX PLANES
  for(int k=lz+2;k<=uz-2;k++)
          for(int i=lx+2;i<=ux-2;i++)
          {
            Phase_field_solver(i,ly+1,k,pfm_dt);
            if(uy-ly>2)Phase_field_solver(i,uy-1,k,pfm_dt);
          }

     cycleArrays(index_grain1,index_grain2);
     cycleArrays(value_grain1,value_grain2);  

 // if(flag_diff)return;
  mpi_config.Start_nonblock_phase(value_phase);
  
  
  for(n=1;n<=num_iter_diff_per_pfm;n++)
  {
    mpi_config.Start_nonblock_conc(total_conc);
      
  #pragma omp parallel for collapse(3)
      // calculate on the core points
    for(int k=lz+2;k<=uz-2;k++)
       for(int j=ly+2;j<=uy-2;j++)
          for(int i=lx+2;i<=ux-2;i++)
          {
              
              Diffusion_solver(i,j,k,new_diff_dt);
          }
    if(n==1)  mpi_config.Wait_nonblock_request(PHASE);
    mpi_config.Wait_nonblock_request(CONC);
     
      //  partition total_conc to value_conc1
    //XY PLANES
    for(int k=lz;k<=uz;k+=(uz-lz))
      for(int j=ly+1;j<=uy-1;j++)
        for(int i=lx+1;i<=ux-1;i++){
          Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
         
          Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);        
    }    
    //YZ PLANES
    for(int k=lz+1;k<=uz-1;k++)
      for(int j=ly+1;j<=uy-1;j++)
        for(int i=lx;i<=ux;i+=(ux-lx)){
          Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
         
          Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);    
    } 
    //ZX PLANES
    for(int k=lz+1;k<=uz-1;k++)
      for(int j=ly;j<=uy;j+=(uy-ly))
        for(int i=lx+1;i<=ux-1;i++){
          Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
         
          Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);    
    }
    // SOLVE DIFFUSION ON SURFACE POINTS
      //XY PLANES    
    for(int j=ly+1;j<=uy-1;j++)
          for(int i=lx+1;i<=ux-1;i++)
          { 
            Diffusion_solver(i,j,lz+1,new_diff_dt);
            if(uz-lz>2)Diffusion_solver(i,j,uz-1,new_diff_dt);
            
          }
      //YZ PLANES   
    for(int k=lz+2;k<=uz-2;k++)
          for(int j=ly+1;j<=uy-1;j++)
          {
            Diffusion_solver(lx+1,j,k,new_diff_dt);
            if(ux-lx>2)  Diffusion_solver(ux-1,j,k,new_diff_dt);  
          }
      //ZX PLANES
    for(int k=lz+2;k<=uz-2;k++)
          for(int i=lx+2;i<=ux-2;i++)
          {
            Diffusion_solver(i,ly+1,k,new_diff_dt);
            if(uy-ly>2)Diffusion_solver(i,uy-1,k,new_diff_dt);
          }
          
    cycleArrays(value_conc1,value_conc2);
     
  }
  return;
}

void phase_field::Update_parameter()
{
    temperature=Temperature_function(time);
    mobility=mobility0*exp(-mobility_Q/(R*temperature));
    C_diffusivity=C_diffusivity0*exp(-C_diffusivity_Q/(R*temperature));
    pfm_dt=dx*dx/dimensions*0.2/(blitz::max(mobility*intf_energy));
    diff_dt=dx*dx/blitz::max(C_diffusivity)/dimensions*0.3;
    thermodynamic.Get_data(temperature);
    time+=pfm_dt;
}

  