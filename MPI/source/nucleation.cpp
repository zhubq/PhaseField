#include "nucleation.h"
#include <cassert>
void nucleation::Nucleate1in2(phase_field & F,const int & num)
  {
    ranlib::Uniform<float> random_number;
    Array<float,3> flag_nucleation(Range(F.lx+1,F.ux-1),Range(F.ly+1,F.uy-1),Range(F.lz+1,F.uz-1),fortranArray);
    flag_nucleation=0.;
    
    int lx=flag_nucleation.lbound(0);
    int ux=flag_nucleation.ubound(0);
    int ly=flag_nucleation.lbound(1);
    int uy=flag_nucleation.ubound(1);
    int lz=flag_nucleation.lbound(2);
    int uz=flag_nucleation.ubound(2);
    int cart_rank=F.mpi_config.Cart_rank();
    MPI_Comm cart_comm=F.mpi_config.Cart_comm();
    random_number.seed(num+cart_rank);
    for(int k=lz;k<=uz;k++)
      for(int j=ly;j<=uy;j++)
        for(int i=lx;i<=ux;i++)
        {
          flag_nucleation(i,j,k)=float(count(F.index_grain1(Range::all(),i,j,k)>0))-float(F.value_conc1(2,i,j,k))*(1.-0.001*random_number.random());
	  
        }
//         assert(any(flag_nucleation>1));
    
   // MPI_allreduce(flag_nucleation.data(),maxLoc,MPI_FLOAT, MPI_MAXLOC,F.mpi_config.Cart_comm());
   // TinyVector<int,3> Extent=flag_nucleation.extent();
   struct{
       float value;
       int rank;
   } local_maxValue,global_maxValue;
   //TinyVector<int,3> maxLocation;
    for(int n=1;n<=num;n++)
    {
      local_maxValue.value=max(flag_nucleation);
      local_maxValue.rank=cart_rank;
      MPI_Allreduce(&local_maxValue,&global_maxValue,1,MPI_FLOAT_INT,MPI_MAXLOC,cart_comm);
      TinyVector<int,3> maxLocation=0;
      if(cart_rank==global_maxValue.rank)
      //{
          maxLocation=maxIndex(flag_nucleation); 
      //    cout<<maxLocation;
          
    //}
      //int a[3];
     // for(int i=0;i<3;i++) a[i]=maxLocation(i);
      MPI_Bcast(maxLocation.data(),3,MPI_INT,global_maxValue.rank,cart_comm);
      //for(int i=0;i<3;i++) maxLocation(i)=a[i];
      //if(cart_rank==0)cout<<"nuclei "<<maxLocation<<endl;
      float radius=1*F.i_intf;
      for(int k=lz;k<=uz;k++)
        for(int j=ly;j<=uy;j++)
          for(int i=lx;i<=ux;i++)
          {
            float distance(TinyVector<int,3> L1, TinyVector<int,3>L2,TinyVector<int,3>Extent);
            float Dist=distance(maxLocation,TinyVector<int,3>(i,j,k),F.domain);
            if(Dist<radius)
            {
              F.index_grain1(Range::all(),i,j,k)=0;
              F.index_grain1(1,i,j,k)=n+100000;
              F.value_grain1(Range::all(),i,j,k)=0.;
              F.value_grain1(1,i,j,k)=1.;
              F.value_phase(Range::all(),i,j,k)=0.;
              F.value_phase(1,i,j,k)=1.;
              //F.value_conc1(1,i,j,k)=F.total_conc(i,j,k);
              flag_nucleation(i,j,k)=0.;
            }else if(Dist<2*radius)
            {
                flag_nucleation(i,j,k)=0.;
            }
            
          }
          
    }
  }
  
float distance(TinyVector<int,3> L1, TinyVector<int,3>L2,TinyVector<int,3>Extent)
{
  TinyVector<float,3>SqrD1=sqr(L1-L2),SqrD2=sqr(Extent-abs(L1-L2));
  float dist=0.;
  for(int i=0;i<3;i++)
  {
    dist+=std::min(SqrD1(i),SqrD2(i));
  }
  return sqrt(dist);
}