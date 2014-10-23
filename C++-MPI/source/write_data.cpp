#include "phase_field.h"
#include <fstream>
#include <iomanip>
void phase_field::Write_data(bool forced )
{
    static double last_temperature=0.;
    static int micro_count=0;
    //int myrank;
    //MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if ( forced||(abs(time-timestep_output*floor(time/timestep_output+0.5))<0.5*pfm_dt)||
          abs(temperature-last_temperature)>10.)
    {        
        Array<float,1> local_vol_fraction(max_num_phase,fortranArray),global_vol_fraction(max_num_phase,fortranArray);
        float domain_vol=ix*iy*iz;
        for(int i=1;i<=max_num_phase;i++){
            local_vol_fraction(i)=sum( value_phase(i,Range(lx+1,ux-1),Range(ly+1,uy-1),Range(lz+1,uz-1)) )/domain_vol;
        
        MPI_Reduce(local_vol_fraction.data(),global_vol_fraction.data(),max_num_phase,MPI_FLOAT,MPI_SUM,0,mpi_config.Cart_comm());
        }
      if(mpi_config.Cart_rank()==0)
      {
        kinetics<<setw(8)<<setprecision(4)<<time<<" "<<temperature<<" ";
        
        for(int i=1;i<=max_num_phase;i++)
             kinetics<<setprecision(4)<<global_vol_fraction(i)<<" ";
        //kinetics<<sum( where(index_grain1<max_num_grains_per_phase, value_grain1, 0.) )/domain_vol;//recrystalized fraction
        kinetics<<endl;
        
        cout<<"time:"<<time<<" "<<"dt: "<<pfm_dt<<" "<<diff_dt<<endl;
      }
      
      if (is_micro_out){
            mpi_config.Start_nonblock_index(index_grain1);
            mpi_config.Start_nonblock_conc(total_conc);
            ostringstream filename;
            //filename<<"field";
            filename.width(7);
            filename.fill('0');
            filename
            //<<int(time*100)<<"_"
	    <<micro_count<<"_"
	    <<mpi_config.Cart_rank()<<".vti";
            ofstream microstructure( (directory+'/'+filename.str()).c_str());
            mpi_config.Wait_nonblock_request(INDEX);
            mpi_config.Wait_nonblock_request(CONC);
            Vti_write(microstructure);
            
            microstructure.close();	  
            
	    if(mpi_config.Cart_rank()==0)
	    {
// #ifdef _DEBUG
//    cout<<"write pvti"<<endl;
// #endif
              ostringstream pfilename;
              //filename<<"field";
              pfilename.width(7);
              pfilename.fill('0');
              pfilename
	      <<micro_count
	      <<".pvti";
              ofstream microstructure( (directory+'/'+pfilename.str()).c_str());
              
              Pvti_write(microstructure,micro_count);
            
              microstructure.close();
	    }
	    micro_count++;
	    
        }
        
        last_temperature=temperature;
    }
    
    
    return;
}
void phase_field::Pvti_write(ofstream & file, int micro_count)
{
      file<<"<?xml version=\"1.0\"?>"<<'\n'
        <<"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<"\n"
        <<"<PImageData WholeExtent=\""<<1<<' '<<ix<<' '
                                   <<1<<' '<<iy<<' '
                                   <<1<<' '<<iz<<'"'<<" GhostLevel=\"1\""
                    <<" Origin="<<"\"0 0 0\" Spacing=\""<<dx<<' '<<dx<<' '<<dx<<"\">"<<endl;
   
     file<<"<PPointData Scalars=\"energy\">"<<endl;
     file<<"<PDataArray type=\"Int32\" Name=\"korn\" format=\"ascii\"/>"<<endl;
     file<<"<PDataArray type=\"Float32\" Name=\"conc\" format=\"ascii\"/>"<<endl;
     file<<"</PPointData>"<<endl;
     int commsize;
     MPI_Comm_size(MPI_COMM_WORLD,&commsize);
     for(int i=0;i<commsize;i++)
     {
            ostringstream filename;
            //filename<<"field";
            filename.width(7);
            filename.fill('0');
            filename
            //<<int(time*100)<<"_"
	    <<micro_count<<"_"
	    <<i<<".vti";
            //ofstream microstructure( (directory+'/'+filename.str()).c_str());
       file<<"<Piece Extent=\""<<max(1,mpi_config.Subdomain_size(i)(0))<<" "<<max(1,mpi_config.Subdomain_size(i)(1)-1)<<" "<<max(1,mpi_config.Subdomain_size(i)(2))
           <<" "<<max(1,mpi_config.Subdomain_size(i)(3)-1)<<" "<<max(1,mpi_config.Subdomain_size(i)(4))<<" "<<max(1,mpi_config.Subdomain_size(i)(5)-1)
	   <<"\" Source=\""<<filename.str()<<"\"/>"<<endl; 
     }
     file<<"</PImageData>"<<endl;
     file<<"</VTKFile>"<<endl;
}
void phase_field::Vti_write(ofstream & file)
{
    int ghost=1;
    int subLx=max(lx,ghost),subLy=max(ly,ghost),subLz=max(lz,ghost);
    int subUx=max(ux-ghost,ghost),subUy=max(uy-ghost,ghost),subUz=max(uz-ghost,ghost);
    file<<"<?xml version=\"1.0\"?>"<<'\n'
        <<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<"\n"
        <<"<ImageData WholeExtent=\""<<subLx<<' '<<subUx<<' '
                                   <<subLy<<' '<<subUy<<' '
                                   <<subLz<<' '<<subUz<<'"'
                    <<" Origin="<<"\"0 0 0\" Spacing=\""<<dx<<' '<<dx<<' '<<dx<<"\">"<<endl;
     file<<"<Piece Extent=\""<<subLx<<' '<<subUx<<' '
                                   <<subLy<<' '<<subUy<<' '
                                   <<subLz<<' '<<subUz<<"\">"<<endl;
     file<<"<PointData Scalars=\"energy\">"<<endl;
     file<<"<DataArray type=\"Int32\" Name=\"korn\" format=\"ascii\">"<<endl;
     for(int k=subLz;k<=subUz;k++)
         for(int j=subLy;j<=subUy;j++)
             for(int i=subLx;i<=subUx;i++)
//                  cout<<i<<" "<<j<<" "<<k<<endl;
                 if(index_grain1(2,i,j,k)>0)
                     file<<0<<" ";
                 else
                     file<<Phase_of_grain(index_grain1(1,i,j,k))<<" ";
                    // file<<(index_grain1(1,i,j,k))<<" ";
                 
     
     file<<"\n</DataArray>"<<endl;
     file<<"<DataArray type=\"Float32\" Name=\"conc\" format=\"ascii\">"<<endl;
     for(int k=subLz;k<=subUz;k++)
         for(int j=subLy;j<=subUy;j++)
             for(int i=subLx;i<=subUx;i++){
                 //file<<setprecision(3)<<value_phase(i,j,k)<<" ";
                 file<<setprecision(3)<<total_conc(i,j,k)<<" ";
                 //cout<<" "<<i<<" "<<j<<" "<<k<<" "<<total_conc(i,j,k)<<endl;
             }
     file<<"</DataArray>"<<endl
         <<"</PointData>"<<endl
         <<"</Piece>"<<endl
         <<"</ImageData>"<<endl
         <<"</VTKFile>"<<endl;
     
}
  