#include "MPI_Class.h"



using namespace std;
mpi_configuration::mpi_configuration()//:subdomain_info(fortranArray)
  {
    int flag;
    MPI_Initialized(&flag);
    if(not flag){
        MPI_Init(NULL,NULL);
    }
    char name[80];int namelen;
    MPI_Get_processor_name(name,&namelen);
   
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    subdomain_info.resize(comm_size);
    if(myrank==0){
         cout<<"The Machine Name is "<<name<<endl;
         cout<<"Total Number of Processes is "<<comm_size<<endl;
         
    }
   // F=pfm;
  }
  
void mpi_configuration::Commit_cart_datatype(const int max_num_coexist,const int max_num_phase)
  {
     int extent[3];
     extent[0]=subdomain_info(myrank)(1)-subdomain_info(myrank)(0)+1;
     extent[1]=subdomain_info(myrank)(3)-subdomain_info(myrank)(2)+1;
     extent[2]=subdomain_info(myrank)(5)-subdomain_info(myrank)(4)+1;
     MPI_Datatype Grain_Ivec, Grain_Dvec,Phase_Dvec;
     
     MPI_Type_contiguous(max_num_coexist,MPI_INT,&Grain_Ivec);
     MPI_Type_commit(&Grain_Ivec);
     MPI_Type_contiguous(max_num_coexist,MPI_DOUBLE,&Grain_Dvec);
     MPI_Type_commit(&Grain_Dvec);
     MPI_Type_contiguous(max_num_phase,MPI_DOUBLE,&Phase_Dvec);
     MPI_Type_commit(&Phase_Dvec);
     
     //     Z
//     ^
//     |       /Y
//     |     /
//     |   /
//     | /
//     |----------------->X
  // XY plane: start from [lx+1,ly+1,lz+1],exclude halo cells   
     MPI_Type_vector(extent[1]-2,(extent[0]-2),extent[0],Grain_Ivec,&data_type[0][0]);
     MPI_Type_commit(&data_type[0][0]);
     MPI_Type_vector(extent[1]-2,(extent[0]-2),extent[0],Grain_Dvec,&data_type[1][0]);
     MPI_Type_commit(&data_type[1][0]);
     MPI_Type_vector(extent[1]-2,(extent[0]-2),extent[0],MPI_DOUBLE,&data_type[2][0]);
     MPI_Type_commit(&data_type[2][0]);
     MPI_Type_vector(extent[1]-2,(extent[0]-2),extent[0],Phase_Dvec,&data_type[3][0]);
     MPI_Type_commit(&data_type[3][0]);    
     
     for(int i=0;i<4;i++)
         data_type[i][1]=data_type[i][0];
     

//     Z
//     ^
//     |       /Y
//     |     /
//     |   /
//     | /
//     |----------------->X     
   //YZ plane: start from [lx,ly,lz], include halo cells  
     MPI_Type_vector(extent[1]*extent[2],1,extent[0],Grain_Ivec,&data_type[0][2]);
     MPI_Type_commit(&data_type[0][2]);
     MPI_Type_vector(extent[1]*extent[2],1,extent[0],Grain_Dvec,&data_type[1][2]);
     MPI_Type_commit(&data_type[1][2]);
     MPI_Type_vector(extent[1]*extent[2],1,extent[0],MPI_DOUBLE,&data_type[2][2]);
     MPI_Type_commit(&data_type[2][2]);
     MPI_Type_vector(extent[1]*extent[2],1,extent[0],Phase_Dvec,&data_type[3][2]);
     MPI_Type_commit(&data_type[3][2]);
     
     
     for(int i=0;i<4;i++)
         data_type[i][3]=data_type[i][2];
          
     
//     Z
//     ^
//     |       /Y
//     |     /
//     |   /
//     | /
//     |----------------->X
     // ZX plane: start from [lx+1,ly+1,lz+1], exclude halo cells
     MPI_Type_vector(extent[2]-2,extent[0]-2,extent[0]*extent[1],Grain_Ivec,&data_type[0][4]);
     MPI_Type_commit(&data_type[0][4]);
     MPI_Type_vector(extent[2]-2,extent[0]-2,extent[0]*extent[1],Grain_Dvec,&data_type[1][4]);
     MPI_Type_commit(&data_type[1][4]);
     MPI_Type_vector(extent[2]-2,extent[0]-2,extent[0]*extent[1],MPI_DOUBLE,&data_type[2][4]);
     MPI_Type_commit(&data_type[2][4]);
     MPI_Type_vector(extent[2]-2,extent[0]-2,extent[0]*extent[1],Phase_Dvec,&data_type[3][4]);
     MPI_Type_commit(&data_type[3][4]);
     
     
     for(int i=0;i<4;i++)
         data_type[i][5]=data_type[i][4];
     

     
     return;
  }
  
void mpi_configuration::Configure_cart(const TinyVector<int,3> & domain)
  {
    #ifdef _DEBUG
   //MPI_Barrier(MPI_COMM_WORLD);cout<<"sdf ";
#endif 
    for(int i=0;i<3;i++)
    {
      cart_periods[i]=1;//periodic condition in each dimension;
      if(domain[i]>1)
          cart_dims[i]=0;
      else
          cart_dims[i]=1;
    }
    
    MPI_Dims_create(comm_size,3,cart_dims);
    for(int i=0;i<3;i++)
    {
        if(domain[i]<cart_dims[i])cout<<"Warning: The number of processes is more than the number of grid points in dimension "<<i<<endl;
    }
    #ifdef _DEBUG
  // MPI_Barrier(MPI_COMM_WORLD);cout<<"sdf ";
#endif 
    
    MPI_Cart_create(MPI_COMM_WORLD,3,cart_dims,cart_periods,0,&cart_comm);
    MPI_Cart_get(cart_comm,3,cart_dims,cart_periods,cart_coord);
   
    //int k=0;
   // for(int i=0;i<3;i++)
      //for(int j=-1;j<=1;j+=2)
      {
         MPI_Cart_shift(cart_comm,2,1,&cart_neighbors[0],&cart_neighbors[1]); // neighbors along Z axis, exchange XY plane
         MPI_Cart_shift(cart_comm,0,1,&cart_neighbors[2],&cart_neighbors[3]);// neigbours along X axis, exchange YZ plane
         MPI_Cart_shift(cart_comm,1,1,&cart_neighbors[4],&cart_neighbors[5]);//neigbours along Y axis, exchange ZX plane
         //k+=2;
      }
#ifdef _DEBUG
  // MPI_Barrier(MPI_COMM_WORLD);cout<<"sdf ";
#endif    
    for(int i=0;i<3;i++)
    {
        int quot=domain[i]/cart_dims[i];
        int mod=domain[i]%cart_dims[i];
        if(cart_coord[i]<mod)
        {
            subdomain_info(myrank)(2*i)=(quot+1)*cart_coord[i];
            subdomain_info(myrank)(2*i+1)=subdomain_info(myrank)(2*i)+(quot+1)+1;
        }
        else
        {
            subdomain_info(myrank)(i*2)=(quot+1)*mod+quot*(cart_coord[i]-mod);
            subdomain_info(myrank)(i*2+1)=subdomain_info(myrank)(i*2)+quot+1;
        }
        int sub_info[6];
        for(int i=0;i<6;i++)
            sub_info[i]=subdomain_info(myrank)(i);

        MPI_Allgather(sub_info,6,MPI_INT,&subdomain_info(0),6,MPI_INT,cart_comm);

    }
    return ;
  }

void mpi_configuration::Start_nonblock_index( Array<int,4> & Index)
  {
    const int i=0;
    int lx=Index.lbound(1);
    int ux=Index.ubound(1);
    int ly=Index.lbound(2);
    int uy=Index.ubound(2);
    int lz=Index.lbound(3);
    int uz=Index.ubound(3);
    int disp;
    int stride;
    void *ptr_base;
  
//     Z
//     ^
//     |       /Y
//     |     /
//     |   /
//     | /
//     |----------------->X
    //XY PLANE 
    
    send_buffer_ptr[i][0]=&Index(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][0]=&Index(1,lx+1,ly+1,uz);      

    //XY PLANE 
    send_buffer_ptr[i][1]=&Index(1,lx+1,ly+1,uz-1);
    recv_buffer_ptr[i][1]=&Index(1,lx+1,ly+1,lz); 

    //YZ PLANE 
    send_buffer_ptr[i][2]=&Index(1,lx+1,ly,lz);
    recv_buffer_ptr[i][2]=&Index(1,ux,ly,lz);

    //YZ PLANE 
    send_buffer_ptr[i][3]=&Index(1,ux-1,ly,lz);
    recv_buffer_ptr[i][3]=&Index(1,lx,ly,lz);

    //ZX PLANE 
    send_buffer_ptr[i][4]=&Index(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][4]=&Index(1,lx+1,uy,lz+1);

    //ZX PLANE 
    send_buffer_ptr[i][5]=&Index(1,lx+1,uy-1,lz+1);
    recv_buffer_ptr[0][5]=&Index(1,lx+1,ly,lz+1);
    

    for(int j=0;j<6;j++)//Index,value_grain1,total_conc
    {   
        nonblock_send_tag[i][j]=myrank;
        MPI_Isend(send_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_send_tag[i][j],cart_comm,&nonblock_send_request[i][j]);
        
        nonblock_recv_tag[i][j]=cart_neighbors[j];
        MPI_Irecv(recv_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_recv_tag[i][j],cart_comm,&nonblock_recv_request[i][j]);

    }

  }
 
void mpi_configuration::Start_nonblock_value( Array<double,4> & Value)
  {
    const int i=1;
    int lx=Value.lbound(1);
    int ux=Value.ubound(1);
    int ly=Value.lbound(2);
    int uy=Value.ubound(2);
    int lz=Value.lbound(3);
    int uz=Value.ubound(3);
  
//     Y
//     ^
//     |          /Z
//     |       /
//     |    /
//     | /
//     |----------------->X
    //XY PLANE FRONT
    send_buffer_ptr[i][0]=&Value(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][0]=&Value(1,lx+1,ly+1,uz);       

    //XY PLANE BACK
    send_buffer_ptr[i][1]=&Value(1,lx+1,ly+1,uz-1);
    recv_buffer_ptr[i][1]=&Value(1,lx+1,ly+1,lz);

    //YZ PLANE LEFT
    send_buffer_ptr[i][2]=&Value(1,lx+1,ly,lz);
    recv_buffer_ptr[i][2]=&Value(1,ux,ly,lz);

    //YZ PLANE RIGHT
    send_buffer_ptr[i][3]=&Value(1,ux-1,ly,lz);
    recv_buffer_ptr[i][3]=&Value(1,lx,ly,lz);

    //ZX PLANE BOTTOM
    send_buffer_ptr[i][4]=&Value(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][4]=&Value(1,lx+1,uy,lz+1);

    //ZX PLANE TOP
    send_buffer_ptr[i][5]=&Value(1,lx+1,uy-1,lz+1);
    recv_buffer_ptr[i][5]=&Value(1,lx+1,ly,lz+1);

    
    //for(int i=0;i<2;i++)// for 6 faces of each subdomain
      for(int j=0;j<6;j++)//Index,value_grain1,total_conc
    {   
        nonblock_send_tag[i][j]=myrank;
        MPI_Isend(send_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_send_tag[i][j],cart_comm,&nonblock_send_request[i][j]);
        
        nonblock_recv_tag[i][j]=cart_neighbors[j];
        MPI_Irecv(recv_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_recv_tag[i][j],cart_comm,&nonblock_recv_request[i][j]);

    }

  }
 
void mpi_configuration::Start_nonblock_phase( Array<double,4> & Phase)
  {
    const int i=3;
    int lx=Phase.lbound(1);
    int ux=Phase.ubound(1);
    int ly=Phase.lbound(2);
    int uy=Phase.ubound(2);
    int lz=Phase.lbound(3);
    int uz=Phase.ubound(3);
  
//     Y
//     ^
//     |          /Z
//     |       /
//     |    /
//     | /
//     |----------------->X
    //XY PLANE FRONT
    send_buffer_ptr[i][0]=&Phase(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][0]=&Phase(1,lx+1,ly+1,uz);     

    //XY PLANE BACK
    send_buffer_ptr[i][1]=&Phase(1,lx+1,ly+1,uz-1);
    recv_buffer_ptr[i][1]=&Phase(1,lx+1,ly+1,lz);  

    //YZ PLANE LEFT
    send_buffer_ptr[i][2]=&Phase(1,lx+1,ly,lz);
    recv_buffer_ptr[i][2]=&Phase(1,ux,ly,lz);

    //YZ PLANE RIGHT
    send_buffer_ptr[i][3]=&Phase(1,ux-1,ly,lz);
    recv_buffer_ptr[i][3]=&Phase(1,lx,ly,lz);

    //ZX PLANE BOTTOM
    send_buffer_ptr[i][4]=&Phase(1,lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][4]=&Phase(1,lx+1,uy,lz+1);

    //ZX PLANE TOP
    send_buffer_ptr[i][5]=&Phase(1,lx+1,uy-1,lz+1);
    recv_buffer_ptr[i][5]=&Phase(1,lx+1,ly,lz+1);
    
    for(int j=0;j<6;j++)//Index,value_grain1,total_conc
    {   
        nonblock_send_tag[i][j]=myrank;
        MPI_Isend(send_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_send_tag[i][j],cart_comm,&nonblock_send_request[i][j]);
        
        nonblock_recv_tag[i][j]=cart_neighbors[j];
        MPI_Irecv(recv_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_recv_tag[i][j],cart_comm,&nonblock_recv_request[i][j]);

    }

  }

void mpi_configuration::Start_nonblock_conc( Array<double,3> & Conc)
  {
    
    const int i=2;
    int lx=Conc.lbound(0);
    int ux=Conc.ubound(0);
    int ly=Conc.lbound(1);
    int uy=Conc.ubound(1);
    int lz=Conc.lbound(2);
    int uz=Conc.ubound(2);
//     Y
//     ^
//     |          /Z
//     |       /
//     |    /
//     | /
//     |----------------->X
    //XY PLANE FRONT
    send_buffer_ptr[i][0]=&Conc(lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][0]=&Conc(lx+1,ly+1,uz);       

    //XY PLANE BACK
    send_buffer_ptr[i][1]=&Conc(lx+1,ly+1,uz-1);
    recv_buffer_ptr[i][1]=&Conc(lx+1,ly+1,lz);

    //YZ PLANE LEFT
    send_buffer_ptr[i][2]=&Conc(lx+1,ly,lz);
    recv_buffer_ptr[i][2]=&Conc(ux,ly,lz);

    //YZ PLANE RIGHT
    send_buffer_ptr[i][3]=&Conc(ux-1,ly,lz);
    recv_buffer_ptr[i][3]=&Conc(lx,ly,lz);

    //ZX PLANE BOTTOM
    send_buffer_ptr[i][4]=&Conc(lx+1,ly+1,lz+1);
    recv_buffer_ptr[i][4]=&Conc(lx+1,uy,lz+1);

    //ZX PLANE TOP
    send_buffer_ptr[i][5]=&Conc(lx+1,uy-1,lz+1);
    recv_buffer_ptr[i][5]=&Conc(lx+1,ly,lz+1);
    
    for(int j=0;j<6;j++)//Index,value_grain1,total_conc
    {   
        nonblock_send_tag[i][j]=myrank;
        MPI_Isend(send_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_send_tag[i][j],cart_comm,&nonblock_send_request[i][j]);
        
        nonblock_recv_tag[i][j]=cart_neighbors[j];
        MPI_Irecv(recv_buffer_ptr[i][j],1,data_type[i][j],cart_neighbors[j],nonblock_recv_tag[i][j],cart_comm,&nonblock_recv_request[i][j]);

    }


  }
  
  
 void mpi_configuration::Free_nonblock_request()
 {
        for(int i=0;i<4;i++)
          for(int j=0;j<6;j++)
        {
            MPI_Request_free(&nonblock_recv_request[i][j]);
            MPI_Request_free(&nonblock_send_request[i][j]);
        }
 }
  
 void mpi_configuration::Wait_nonblock_request(int i)
 {   
   // for(int j=0;j<6;i++)
    {
       MPI_Waitall(6,nonblock_recv_request[i],nonblock_recv_status[i]);
       MPI_Waitall(6,nonblock_send_request[i],nonblock_send_status[i]);
    }
 }
 
   

  