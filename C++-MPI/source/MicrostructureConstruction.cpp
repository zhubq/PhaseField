#include "phase_field.h"
#include "random/uniform.h"
#include "MPI_Class.h"
void phase_field::Construct(bool c,char a[])
{
    if (c){
      cout<<"using Voronoi tessellation to construct the structure."<<endl;  
      Array<int, 3> grain_index(Range(lx,ux),Range(ly,uy),Range(lz,uz),fortranArray);//used for tessellation
    
      Tessellation(grain_index,num_grain[1],reduction);//use Voronoi tessellation to generate a polycrystalline microstructure

      for(int k=lz;k<=uz;k++)
          for(int j=ly;j<=uy;j++)
              for(int i=lx;i<=ux;i++)		
                  index_grain1(1,i,j,k)=grain_index(i,j,k)+200000; // get a austenitic structure
    }else{
            cout<<"read initial microstructure from file "<<a<<endl;
            ifstream file(a);
	    for(int k=lz;k<=uz;k++)
	      for(int j=ly;j<=uy;j++)
		for(int i=lx;i<=ux;i++)
	              file>>index_grain1(1,i,j,k);
            file.close();            
    }
    value_grain1(1,Range::all(),Range::all(),Range::all())=1.0; // assign volume fraction of each grain: 1.0 at the beginning without diffuse interface

#pragma omp parallel for collapse(3)
    for(int k=lz;k<=uz;k++)
      for(int j=ly;j<=uy;j++)
        for (int i=lx;i<=ux;i++)
        {  
            
            value_phase(Phase_of_grain(index_grain1(1,i,j,k)),i,j,k)=1.0;
            
        } // assign volume fraction of each phase
    for(int k=lz;k<=uz;k++)
        for(int j=ly;j<=uy;j++)
            for(int i=lx;i<=ux;i++)     
                for(int n=0;n<max_num_phase;n++)
                    value_conc1(n+1,i,j,k)=init_conc[n];// initialize phase concentrations
    // cout<<value_conc1(1,lx,ly,lz)<<endl;
    /* GENERATE DIFFUSE INTERFACES*/
    if(num_init_interation>0)cout<<"start initial iteration to generate diffuse interfaces"<<endl;
 
    for(int n=1;n<=num_init_interation;n++) {
      // send-recv halo cells  
      mpi_config.Start_nonblock_index(index_grain1);
      mpi_config.Start_nonblock_value(value_grain1);
      
  #pragma omp parallel for collapse(3)
      // calculate on the core points
      for(int k=lz+2;k<=uz-2;k++)
        for(int j=ly+2;j<=uy-2;j++)
          for(int i=lx+2;i<=ux-2;i++)
          {
              
              Iso_growth(i,j,k);
          }
      mpi_config.Wait_nonblock_request(0);
      mpi_config.Wait_nonblock_request(1);

      // calculate surface points
      //XY PLANES
      for(int j=ly+1;j<=uy-1;j++)
          for(int i=lx+1;i<=ux-1;i++)
          {   //cout<<"clc ";
              Iso_growth(i,j,lz+1);
              if(uz-lz>2)Iso_growth(i,j,uz-1);
          }

      //YZ PLANES
      for(int k=lz+2;k<=uz-2;k++)
          for(int j=ly+1;j<=uy-1;j++)
          {
            Iso_growth(lx+1,j,k);
            if(ux-lx>2)Iso_growth(ux-1,j,k);
          }
      //ZX PLANES
      for(int k=lz+2;k<=uz-2;k++)
          for(int i=lx+2;i<=ux-2;i++)
          {
            Iso_growth(i,ly+1,k);
            if(uy-ly>2)Iso_growth(i,uy-1,k);
          }

     cycleArrays(index_grain1,index_grain2);
     cycleArrays(value_grain1,value_grain2);
     
    }
   
  Update_parameter();

  mpi_config.Start_nonblock_phase(value_phase);
  mpi_config.Wait_nonblock_request(PHASE);
 
#pragma omp for
    for(int k=lz;k<=uz;k++)
      for(int j=ly;j<=uy;j++)
        for(int i=lx;i<=ux;i++){
	  total_conc(i,j,k)=dot( value_phase(Range::all(),i,j,k),value_conc1(Range::all(),i,j,k) );//cout<<" "<<i<<" "<<j<<" "<<k<<" "<<value_phase(Range::all(),i,j,k)<<" "<<total_conc(i,j,k)<<endl;
          Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
         
          Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);
    //exit(0);
    
    }
    //cout<<total_conc(lx+1,ly+1,lz+1)<<endl;
//     for(int k=lz;k<=uz;k++)
//       for(int j=ly;j<=uy;j++)
//         for(int i=lx;i<=ux;i++)
//             cout<<" "<<i<<" "<<j<<" "<<k<<" "<<total_conc(i,j,k)<<endl;
    
    if(num_init_interation>0)cout<<"******Diffuse Interface Generated******"<<endl;//cout<<pi*pi/(2.*eta*eta)*(dx*dx)<<endl;
}


    
void   phase_field::Tessellation(Array<int,3> & grain_index,const int num, const double reduction)
{
    //coordinates of grain center for tessellation.
    Array<double,2> point_xyz(3,num,fortranArray);
    //generate normalized [0,0.1] random coordinates by Master process
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank==0)
    {
      ranlib::UniformClosed<double> uniform_numbers;
      uniform_numbers.seed(num);
  
      for(Array<double,2>::iterator i=point_xyz.begin(); i!=point_xyz.end();i++)
      {
         *i=uniform_numbers.random();
      }
    }
    // Broadcast the coordinates to all processes from Process 0[master]
    MPI_Bcast(point_xyz.data(),3*num,MPI_DOUBLE,0,mpi_config.Cart_comm());
    //Scale coordinates to the domain size
    point_xyz(1,Range::all())*=(ix-1);
    point_xyz(1,Range::all())+=1;
    
    point_xyz(2,Range::all())*=(iy-1);
    point_xyz(2,Range::all())+=1;
    
    point_xyz(3,Range::all())*=(iz-1);
    point_xyz(3,Range::all())+=1; 
    
    
    //TinyVector<int,3> lbound=grain_index.lbound(), ubound=grain_index.ubound(),shape=grain_index.shape();
#pragma omp parallel for //shared(grain_index,num)
    for(int k=lz;k<=uz;k++)
      for(int j=ly;j<=uy;j++)
	for(int i=lx;i<=ux;i++){
	  int index=0;
	  double temp1=dot(domain/(1.-reduction),domain/(1.-reduction));
	  double temp2;
	  for(int n=1;n<=num;n++)
	  {
	    temp2=pow2( ( abs(point_xyz(1,n)-i)<0.5*ix?abs(point_xyz(1,n)-i):ix-abs(point_xyz(1,n)-i) )/(1.-reduction) )
	                +pow2( ( abs(point_xyz(2,n)-j)<0.5*iy?abs(point_xyz(2,n)-j):iy-abs(point_xyz(2,n)-j) )*(1.-reduction) )
			+pow2( abs(point_xyz(3,n)-k)<0.5*iz?abs(point_xyz(3,n)-k):iz-abs(point_xyz(3,n)-k) );
	    if(temp1>temp2){
	      temp1=temp2;
	      index=n;
	    }
	  }
	  grain_index(i,j,k)=index;
	}
	
	  
}

static void Select_bands(Array<int,3> & grain_index,int num,double fraction,int spacing)
{
  int num_pixel=grain_index.size()*fraction;
  Array<int,1> grain_volume(num,fortranArray);
  
  for(int i=1;i<=num;i++)
    grain_volume(i)=count(grain_index==i);
  
  Array<int,1> Is_InBands(num,fortranArray);
#pragma omp for
  for (int k=grain_index.lbound(2);k<=grain_index.ubound(2);k++)
    for(int j=grain_index.lbound(1)+1;j<grain_index.ubound(1);j+=spacing)
      for(int i=grain_index(0);i<=grain_index.ubound(0);i++)
      {
	++Is_InBands(grain_index(i-1,j,k));
	++Is_InBands(grain_index(i,j,k));
        ++Is_InBands(grain_index(i+1,j,k));
      }
      
} 

// void phase_field::Iso_growth(const int i,const int j,const int k)
// {//cout<<"start"<<endl;
//   const static int max_num_nonzero_stencil=10;
//   Array<int,1> index_grain_nonzero_stencil(max_num_nonzero_stencil,fortranArray);//assume there are at most 10 non-zero phase-fields on all stencil-points,7-stencil pattern is used in 3D.
//   index_grain_nonzero_stencil=0;
//   Array<int,1> index_phase_nonzero_stencil(max_num_nonzero_stencil,fortranArray);
//   Array<double,2> value_grain_nonzero_stencil(max_num_nonzero_stencil,2,fortranArray);
//   Array<double,1> del_value(max_num_nonzero_stencil,fortranArray);
//   Array<double,1> local_conc(value_conc1(Range::all(),i,j,k));
//   Array<int,2> neighbors(3,6,fortranArray); 
//   
//   int num_nonzero_local=0,num_nonzero_stencil=0;//number of nonzeros in index_grain_nonzero_stencil.
//  
//   index_grain_nonzero_stencil=0;
//   value_grain_nonzero_stencil=0.;
//   del_value=0.;
//   
//   neighbors=
//              Get_neighbor(i,ix,-1),j,k,
//              Get_neighbor(i,ix,1 ),j,k,
//              i,Get_neighbor(j,iy,-1),k,
//              i,Get_neighbor(j,iy,1 ),k,
//              i,j,Get_neighbor(k,iz,-1),
//              i,j,Get_neighbor(k,iz,1 );
//  
//   for(int ii=1;ii<=max_num_coexist;ii++){
//     if(value_grain1(ii,i,j,k)>cutoff){
//       num_nonzero_stencil++;
//       index_grain_nonzero_stencil(ii)=index_grain1(ii,i,j,k);
//       value_grain_nonzero_stencil(ii,1)=value_grain1(ii,i,j,k);
//     }      
//   }
//   
//   //num_nonzero_stencil=num_nonzero_local;
// 
//   for(int ii=1;ii<=2*dimensions;ii++){
//     for(int jj=1;jj<=max_num_coexist;jj++){
//         int temp_index=index_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
//         if(temp_index==0)break;
//         double temp_value=value_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
// 	int kk;
// 	if(temp_value>cutoff){
//           for(kk=1;kk<=max_num_nonzero_stencil;kk++){
//               if(temp_index==index_grain_nonzero_stencil(kk))
// 	      {value_grain_nonzero_stencil(kk,2)+=temp_value;break;}
// 	  }
// 	}
// 	if(kk>max_num_nonzero_stencil && num_nonzero_stencil<max_num_nonzero_stencil){
//    
//             num_nonzero_stencil+=1;   
//             index_grain_nonzero_stencil(num_nonzero_stencil)=temp_index;               
//             value_grain_nonzero_stencil(num_nonzero_stencil,2)+=temp_value;
//             
//         }
//     }
//   }
// 
//   if(num_nonzero_stencil==1){
//     index_grain2(Range::all(),i,j,k)=index_grain1(Range::all(),i,j,k);
//     value_grain2(Range::all(),i,j,k)=value_grain1(Range::all(),i,j,k);  
//     return;
//   }
//  
// //---------------------------------------------------------------------------------------------
// //=============================================================================================
// /*======now,three arrays to be used for calculation: 
//  * index_grain_nonzero_stencil,index_phase_nonzero_stencil,value_grain_nonzero_stencil.
//  * 
// */
// 
//   del_value=0.;
//   double factor1=1./(dx*dx),
//          factor2=pi*pi/(2.*eta*eta);
// 
//   double del_value_ij;       
//   for(int ii=1;ii<=num_nonzero_stencil;ii++){
// 
//       for(int jj=1;(jj<=num_nonzero_stencil);jj++){
//         if (ii==jj)continue;
// 
//         
//         del_value_ij=( value_grain_nonzero_stencil(ii,2)*value_grain_nonzero_stencil(jj,1)-
//                        value_grain_nonzero_stencil(jj,2)*value_grain_nonzero_stencil(ii,1) )*factor1
//                       +(value_grain_nonzero_stencil(ii,1)-value_grain_nonzero_stencil(jj,1))*factor2;
// 
//         del_value(ii)+=del_value_ij;
//       }
//       del_value(ii)*=0.5/dimensions*dx*dx;
//       del_value(ii)+=value_grain_nonzero_stencil(ii,1);
// 
//       if(del_value(ii)>1.0)del_value(ii)=1.0;
//           else if(del_value(ii)<cutoff)del_value(ii)=0.;
//       
//   }
// 
//   index_grain2(Range::all(),i,j,k)=0;
//   value_grain2(Range::all(),i,j,k)=0.;
//   value_phase(Range::all(),i,j,k)=0.;
// 
//   Sort_xyz(del_value,index_grain_nonzero_stencil,max_num_nonzero_stencil);
// 
//   double sum_del_value=0.;
// 
//   for(int ii=1;ii<=max_num_coexist;ii++){
//     if (del_value(ii)>cutoff){
//         index_grain2(ii,i,j,k)=index_grain_nonzero_stencil(ii);
//         value_grain2(ii,i,j,k)=del_value(ii);
//         int index=Phase_of_grain( index_grain_nonzero_stencil(ii) );
//         value_phase(index,i,j,k)+=del_value(ii);
//         sum_del_value+=del_value(ii);
//     }
//   }
// 
//   value_grain2(Range::all(),i,j,k)*=1./sum_del_value;
//   value_phase(Range::all(),i,j,k)*=1./sum_del_value;
// 
//   return;
// }

void phase_field::Iso_growth(const int i,const int j,const int k)
{//cout<<"start"<<endl;
  const static int max_num_nonzero_stencil=10;
  Array<int,1> index_grain_nonzero_stencil(max_num_nonzero_stencil,fortranArray);//assume there are at most 10 non-zero phase-fields on all stencil-points,7-stencil pattern is used in 3D.
  index_grain_nonzero_stencil=0;
  Array<int,1> index_phase_nonzero_stencil(max_num_nonzero_stencil,fortranArray);
  Array<double,2> value_grain_nonzero_stencil(max_num_nonzero_stencil,2,fortranArray);
  Array<double,1> del_value(max_num_nonzero_stencil,fortranArray);
  Array<double,1> local_conc(value_conc1(Range::all(),i,j,k));
  Array<int,2> neighbors(3,6,fortranArray); 
  
  int num_nonzero_local=0,num_nonzero_stencil=0;//number of nonzeros in index_grain_nonzero_stencil.
 
  index_grain_nonzero_stencil=0;
  value_grain_nonzero_stencil=0.;
  del_value=0.;
  
  neighbors=
             Get_neighbor(i,ix,-1),j,k,
             Get_neighbor(i,ix,1 ),j,k,
             i,Get_neighbor(j,iy,-1),k,
             i,Get_neighbor(j,iy,1 ),k,
             i,j,Get_neighbor(k,iz,-1),
             i,j,Get_neighbor(k,iz,1 );
 
  for(int ii=1;ii<=max_num_coexist;ii++){
    if(value_grain1(ii,i,j,k)>cutoff){
      num_nonzero_stencil++;
      index_grain_nonzero_stencil(ii)=index_grain1(ii,i,j,k);
      value_grain_nonzero_stencil(ii,1)=value_grain1(ii,i,j,k);
    }      
  }
  
  //num_nonzero_stencil=num_nonzero_local;

  for(int ii=1;ii<=2*dimensions;ii++){
    for(int jj=1;jj<=max_num_coexist;jj++){
        int temp_index=index_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
        if(temp_index==0)break;
        double temp_value=value_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
        int kk;
        if(temp_value>cutoff){
          for(kk=1;kk<=max_num_nonzero_stencil;kk++){
              if(temp_index==index_grain_nonzero_stencil(kk))
              {value_grain_nonzero_stencil(kk,2)+=temp_value;break;}
          }
        }
        if(kk>max_num_nonzero_stencil && num_nonzero_stencil<max_num_nonzero_stencil){
   
            num_nonzero_stencil+=1;   
            index_grain_nonzero_stencil(num_nonzero_stencil)=temp_index;               
            value_grain_nonzero_stencil(num_nonzero_stencil,2)+=temp_value;
            
        }
    }
  }

  if(num_nonzero_stencil==1){
    index_grain2(Range::all(),i,j,k)=index_grain1(Range::all(),i,j,k);
    value_grain2(Range::all(),i,j,k)=value_grain1(Range::all(),i,j,k);  
    return;
  }
 
//---------------------------------------------------------------------------------------------
//=============================================================================================
/*======now,three arrays to be used for calculation: 
 * index_grain_nonzero_stencil,index_phase_nonzero_stencil,value_grain_nonzero_stencil.
 * 
*/

  del_value=0.;
  double factor1=1./(dx*dx),
         factor2=pi*pi/(2.*eta*eta);

  double del_value_ij;       
  for(int ii=1;ii<=num_nonzero_stencil;ii++){

      for(int jj=1;(jj<=num_nonzero_stencil);jj++){
        if (ii==jj)continue;

        
        del_value_ij=( value_grain_nonzero_stencil(ii,2)*value_grain_nonzero_stencil(jj,1)-
                       value_grain_nonzero_stencil(jj,2)*value_grain_nonzero_stencil(ii,1) )*factor1
                      +(value_grain_nonzero_stencil(ii,1)-value_grain_nonzero_stencil(jj,1))*factor2;

        del_value(ii)+=del_value_ij;
      }
      del_value(ii)*=0.1/dimensions*dx*dx;
      del_value(ii)+=value_grain_nonzero_stencil(ii,1);

      if(del_value(ii)>1.0)del_value(ii)=1.0;
          else if(del_value(ii)<cutoff)del_value(ii)=0.;
      
  }

  index_grain2(Range::all(),i,j,k)=0;
  value_grain2(Range::all(),i,j,k)=0.;
  value_phase(Range::all(),i,j,k)=0.;

Sort_xyz(del_value,index_grain_nonzero_stencil,max_num_nonzero_stencil);

  double sum_del_value=0.;

  for(int ii=1;ii<=max_num_coexist;ii++){
    if (del_value(ii)>cutoff){
        index_grain2(ii,i,j,k)=index_grain_nonzero_stencil(ii);
        value_grain2(ii,i,j,k)=del_value(ii);
        int index=Phase_of_grain( index_grain_nonzero_stencil(ii) );
        value_phase(index,i,j,k)+=del_value(ii);
        sum_del_value+=del_value(ii);
    }
  }
#ifdef _NDEBUG
 assert(sum_del_value>0.5)
#endif
  value_grain2(Range::all(),i,j,k)*=1./sum_del_value;
  value_phase(Range::all(),i,j,k)*=1./sum_del_value;

  return;
}