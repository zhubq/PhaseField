#include "phase_field.h"
const static int max_num_nonzero_stencil=10;
// void phase_field::Phase_field_solver(const int i,const int j,const int k,const real Pfm_dt) {
//     Array<int,1> index_grain_nonzero_stencil(max_num_nonzero_stencil,fortranArray);//assume there are at most 10 non-zero phase-fields on all stencil-points,7-stencil pattern is used in 3D.
//     index_grain_nonzero_stencil=0;
//     Array<int,1> index_phase_nonzero_stencil(max_num_nonzero_stencil,fortranArray);
//     index_phase_nonzero_stencil=0;
//     int num_nonzero_local=0,num_nonzero_stencil=0;//number of nonzeros in index_grain_nonzero_stencil.
// // int neighbors[7][3]={0};
// //   int ii,jj,kk,nn;
//     Array<real,2> value_grain_nonzero_stencil(max_num_nonzero_stencil,2,fortranArray);
//     value_grain_nonzero_stencil=0.;
//     //cout<<"call"<<endl;
//     Array<int,2> neighbors(3,6,fortranArray);
//     neighbors=
//         Get_neighbor(i,ix,-1),j,k,
//         Get_neighbor(i,ix,1 ),j,k,
//         i,Get_neighbor(j,iy,-1),k,
//         i,Get_neighbor(j,iy,1 ),k,
//         i,j,Get_neighbor(k,iz,-1),
//         i,j,Get_neighbor(k,iz,1 );
//
//     for(int ii=1; ii<=max_num_coexist; ii++) {
//         if(value_grain1(ii,i,j,k)>cutoff) {
//             num_nonzero_local+=1;
//             index_grain_nonzero_stencil(ii)=index_grain1(ii,i,j,k);
//             index_phase_nonzero_stencil(ii)=Phase_of_grain(index_grain1(ii,i,j,k));
//             value_grain_nonzero_stencil(ii,1)=value_grain1(ii,i,j,k);
//         }
//     }
//
//     num_nonzero_stencil=num_nonzero_local;
//
//     for(int ii=1; ii<=2*dimensions; ii++) {
//         for(int jj=1; jj<=max_num_coexist; jj++) {
//             int temp_index=index_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
//
//             if(temp_index==0)break;
//
//             //cout<<"d";
//             real temp_value=value_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
//             int kk;
//
//             for(kk=1; kk<=max_num_nonzero_stencil; kk++) {
//                 if(temp_index==index_grain_nonzero_stencil(kk)) {
//                     value_grain_nonzero_stencil(kk,2)+=temp_value;
//                     break;
//                 }
//             }
//
//             //cout<<kk<<" "<<max_num_nonzero_stencil<<endl;
//             if(kk>max_num_nonzero_stencil&&num_nonzero_stencil<max_num_nonzero_stencil) {
//                 //cout<<"d";
//                 num_nonzero_stencil+=1;
//                 index_grain_nonzero_stencil(num_nonzero_stencil)=temp_index;
//                 index_phase_nonzero_stencil(num_nonzero_stencil)=Phase_of_grain(temp_index);
//                 value_grain_nonzero_stencil(num_nonzero_stencil,2)+=temp_value;
//             }
//         }
//     }
//
//     if(num_nonzero_stencil==1) {
//         for(int ii=1; ii<=max_num_coexist; ii++) {
//             index_grain2(ii,i,j,k)=index_grain1(ii,i,j,k);
//             value_grain2(ii,i,j,k)=value_grain1(ii,i,j,k);
//         }
//
//         return;
//     }
//
// //---------------------------------------------------------------------------------------------
// //=============================================================================================
//     /*======now,three arrays to be used for calculation:
//      * index_grain_nonzero_stencil,index_phase_nonzero_stencil,value_grain_nonzero_stencil.
//      *
//     */
//     Array<real,1> del_value(max_num_nonzero_stencil,fortranArray);
//     del_value=0.;
//     real factor1=1./(dx*dx),
//          factor2=pi*pi/(2.*eta*eta),
//          factor3=pi/eta;
//
//     for(int ii=1; ii<=num_nonzero_stencil; ii++) {
//         int grainii=index_grain_nonzero_stencil(ii);
//         int phaseii=index_phase_nonzero_stencil(ii);
//
//         for(int jj=1; (jj<=num_nonzero_stencil); jj++) {
//             if (ii==jj)continue;
//
//             int grainjj=index_grain_nonzero_stencil(jj);
//             int phasejj=index_phase_nonzero_stencil(jj);
//
//             // if(grainii<max_num_grains_per_phase&&grainjj<max_num_grains_per_phase)continue;//no interaction between deformed grain.
//             if(abs(phaseii-phasejj)==2)continue; //no interaction between pearlite and ferrite.
//
//             real del_value_ij=0.;
//             del_value_ij=(value_grain_nonzero_stencil(ii,2)*value_grain_nonzero_stencil(jj,1)-
//                           value_grain_nonzero_stencil(jj,2)*value_grain_nonzero_stencil(ii,1))*factor1;
//             del_value_ij+=(value_grain_nonzero_stencil(ii,1)-value_grain_nonzero_stencil(jj,1))*factor2;
//             del_value_ij*=intf_energy(phaseii,phasejj);
// //
//             del_value_ij+=factor3*sqrt( value_grain_nonzero_stencil(ii,1)*value_grain_nonzero_stencil(jj,1) )*
//                           Driving_force(grainii,phaseii,grainjj,phasejj,value_conc1(phaseii,i,j,k),value_conc1(phasejj,i,j,k));
//             del_value_ij*=mobility(phaseii,phasejj);
//             del_value(ii)+=del_value_ij;
//         }
//
//         del_value(ii)*=Pfm_dt;
//         del_value(ii)+=value_grain_nonzero_stencil(ii,1);
//
//         if(del_value(ii)>1.0)del_value(ii)=1.0;
//
//         else if(del_value(ii)<cutoff)del_value(ii)=0.;
//     }
//
//     Sort_xyz(del_value,index_grain_nonzero_stencil,max_num_nonzero_stencil);
//     real sum_del_value=0.;
//
//     for(int ii=1; ii<=max_num_coexist; ii++) {
//         int index=Phase_of_grain( index_grain_nonzero_stencil(ii) );
//
//         if (del_value(ii)>cutoff) {
//             index_grain2(ii,i,j,k)=index_grain_nonzero_stencil(ii);
//             value_grain2(ii,i,j,k)=del_value(ii);
//             value_phase(index,i,j,k)+=del_value(ii);
//             sum_del_value+=del_value(ii);
//
//         } else {
//             index_grain2(ii,i,j,k)=0;
//             value_grain2(ii,i,j,k)=0.;
//             value_phase(index,i,j,k)=0.;
//         }
//     }
//
//     for(int ii=1; ii<=max_num_coexist; ii++)
//         value_grain2(ii,i,j,k)*=1./sum_del_value;
//
//     for(int ii=1; ii<=max_num_phase; ii++)
//         value_phase(ii,i,j,k)*=1./sum_del_value;
//
//     //Partitioning(total_conc(i,j,k),&value_phase(1,i,j,k),&value_conc1(1,i,j,k));
//     Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
//     Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);
//     return;
// }

void phase_field::Phase_field_solver(const int i,const int j,const int k,const double Pfm_dt) {
    const static int max_num_nonzero_stencil=10;
    Array<int,1> index_grain_nonzero_stencil(max_num_nonzero_stencil,fortranArray);//assume there are at most 10 non-zero phase-fields on all stencil-points,7-stencil pattern is used in 3D.
    index_grain_nonzero_stencil=0;
    Array<int,1> index_phase_nonzero_stencil(max_num_nonzero_stencil,fortranArray);
    index_phase_nonzero_stencil=0;
    int num_nonzero_local=0,num_nonzero_stencil=0;//number of nonzeros in index_grain_nonzero_stencil.
// int neighbors[7][3]={0};
//   int ii,jj,kk,nn;
    Array<double,2> value_grain_nonzero_stencil(max_num_nonzero_stencil,2,fortranArray);
    value_grain_nonzero_stencil=0.;
    //cout<<"call"<<endl;
    Array<int,2> neighbors(3,6,fortranArray);
    neighbors=
        Get_neighbor(i,ix,-1),j,k,
        Get_neighbor(i,ix,1 ),j,k,
        i,Get_neighbor(j,iy,-1),k,
        i,Get_neighbor(j,iy,1 ),k,
        i,j,Get_neighbor(k,iz,-1),
        i,j,Get_neighbor(k,iz,1 );

    for(int ii=1; ii<=max_num_coexist; ii++) {
        if(value_grain1(ii,i,j,k)>cutoff) {
            num_nonzero_local+=1;
            index_grain_nonzero_stencil(ii)=index_grain1(ii,i,j,k);
            index_phase_nonzero_stencil(ii)=Phase_of_grain(index_grain1(ii,i,j,k));
            value_grain_nonzero_stencil(ii,1)=value_grain1(ii,i,j,k);
        }
    }

    num_nonzero_stencil=num_nonzero_local;

    for(int ii=1; ii<=2*dimensions; ii++) {
        for(int jj=1; jj<=max_num_coexist; jj++) {
            int temp_index=index_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));

            if(temp_index==0)break;

            //cout<<"d";
            double temp_value=value_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
            int kk;

            for(kk=1; kk<=max_num_nonzero_stencil; kk++) {
                if(temp_index==index_grain_nonzero_stencil(kk)) {
                    value_grain_nonzero_stencil(kk,2)+=temp_value;
                    break;
                }
            }

            //cout<<kk<<" "<<max_num_nonzero_stencil<<endl;
            if(kk>max_num_nonzero_stencil&&num_nonzero_stencil<max_num_nonzero_stencil) {
                //cout<<"d";
                num_nonzero_stencil+=1;
                index_grain_nonzero_stencil(num_nonzero_stencil)=temp_index;
                index_phase_nonzero_stencil(num_nonzero_stencil)=Phase_of_grain(temp_index);
                value_grain_nonzero_stencil(num_nonzero_stencil,2)+=temp_value;
            }
        }
    }

    if(num_nonzero_stencil==1) {
        for(int ii=1; ii<=max_num_coexist; ii++) {
            index_grain2(ii,i,j,k)=index_grain1(ii,i,j,k);
            value_grain2(ii,i,j,k)=value_grain1(ii,i,j,k);
        }

        // index_grain2(Range::all(),i,j,k)=index_grain1(Range::all(),i,j,k);
        //value_grain2(Range::all(),i,j,k)=value_grain1(Range::all(),i,j,k);
        return;
    }

//---------------------------------------------------------------------------------------------
//=============================================================================================
    /*======now,three arrays to be used for calculation:
     * index_grain_nonzero_stencil,index_phase_nonzero_stencil,value_grain_nonzero_stencil.
     *
    */
    Array<double,1> del_value(max_num_nonzero_stencil,fortranArray);
    del_value=0.;
    double factor1=1./(dx*dx),
           factor2=pi*pi/(2.*eta*eta),
           factor3=pi/eta;
    Array<double,1> local_conc(max_num_phase,fortranArray);

    for(int ii=1; ii<=max_num_phase; ii++)
        local_conc(ii)=value_conc1(ii,i,j,k);

    //cout<<"num_non "<<num_nonzero_stencil<<endl;
    for(int ii=1; ii<=num_nonzero_stencil; ii++) {
        int grainii=index_grain_nonzero_stencil(ii);
        int phaseii=index_phase_nonzero_stencil(ii);

        for(int jj=1; (jj<=num_nonzero_stencil); jj++) {
            if (ii==jj)continue;

            int grainjj=index_grain_nonzero_stencil(jj);
            int phasejj=index_phase_nonzero_stencil(jj);

            // if(grainii<max_num_grains_per_phase&&grainjj<max_num_grains_per_phase)continue;//no interaction between deformed grain.
            if(abs(phaseii-phasejj)==2)continue; //no interaction between pearlite and ferrite.

            double del_value_ij=0.;
            del_value_ij=(value_grain_nonzero_stencil(ii,2)*value_grain_nonzero_stencil(jj,1)-
                          value_grain_nonzero_stencil(jj,2)*value_grain_nonzero_stencil(ii,1))*factor1;
            del_value_ij+=(value_grain_nonzero_stencil(ii,1)-value_grain_nonzero_stencil(jj,1))*factor2;
            del_value_ij*=intf_energy(phaseii,phasejj);
//
            del_value_ij+=factor3*sqrt( value_grain_nonzero_stencil(ii,1)*value_grain_nonzero_stencil(jj,1) )*
                          Driving_force(grainii,phaseii,grainjj,phasejj,local_conc(phaseii),local_conc(phasejj));
            del_value_ij*=mobility(phaseii,phasejj);
            del_value(ii)+=del_value_ij;
        }

        del_value(ii)*=Pfm_dt;
        del_value(ii)+=value_grain_nonzero_stencil(ii,1);

        if(del_value(ii)>1.0)del_value(ii)=1.0;

        else if(del_value(ii)<cutoff)del_value(ii)=0.;
    }

    //index_grain2(Range::all(),i,j,k)=0;
    //value_grain2(Range::all(),i,j,k)=0.;
    for(int ii=1; ii<=max_num_phase; ii++)
        value_phase(ii,i,j,k)=0.;

    Sort_xyz(del_value,index_grain_nonzero_stencil,max_num_nonzero_stencil);
// for(int ii=1;ii<=max_num_nonzero_stencil;ii++)
    //   index_phase_nonzero_stencil(ii)=Phase_of_grain(index_grain_nonzero_stencil(ii));
    double sum_del_value=0.;

    for(int ii=1; ii<=max_num_coexist; ii++) {
        if (del_value(ii)>cutoff) {
            index_grain2(ii,i,j,k)=index_grain_nonzero_stencil(ii);
            value_grain2(ii,i,j,k)=del_value(ii);
            int index=Phase_of_grain( index_grain_nonzero_stencil(ii) );
            value_phase(index,i,j,k)+=del_value(ii);
            sum_del_value+=del_value(ii);

        } else {
            index_grain2(ii,i,j,k)=0;
            value_grain2(ii,i,j,k)=0.;
        }
    }

    for(int ii=1; ii<=max_num_coexist; ii++)
        value_grain2(ii,i,j,k)*=1./sum_del_value;

    for(int ii=1; ii<=max_num_phase; ii++)
        value_phase(ii,i,j,k)*=1./sum_del_value;

    //Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),local_conc);
    Partitioning(total_conc(i,j,k),&value_phase(1,i,j,k),&value_conc1(1,i,j,k));
    //value_conc2(Range::all(),i,j,k)//=local_conc;
// cout<<"run"<<endl;
    ;
    return;
}

/*
void phase_field::Diffusion_solver(const int i,const int j,const int k,const double diff_dt)
{
  //Array<int,2> neighbors(3,6,fortranArray);
  int neighbors[6][3];
#pragma ivdep
  for(int ii=0;ii<6;ii++)
  {
    neighbors[ii][0]=i;
    neighbors[ii][1]=j;
    neighbors[ii][2]=k;
  }
  neighbors[0][0]=Get_neighbor(i,ix,-1);
  neighbors[1][0]=Get_neighbor(i,ix,1 );
  neighbors[2][1]=Get_neighbor(j,iy,-1);
  neighbors[3][1]=Get_neighbor(j,iy,1 );
  neighbors[4][2]=Get_neighbor(k,iz,-1);
  neighbors[5][2]=Get_neighbor(k,iz,1 );

  double del_conc=0.;
  for(int ii=0;ii<2*dimensions;ii++)
      for(int jj=1;jj<=max_num_phase;jj++)
         // del_conc+=sqrt( value_phase(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii))*value_phase(jj,i,j,k) )*
         //               ( value_conc1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii))-value_conc1(jj,i,j,k) )*C_diffusivity(jj);

	  del_conc+=sqrt( value_phase(jj,neighbors[ii][0],neighbors[ii][1],neighbors[ii][2])*value_phase(jj,i,j,k) )*
                        ( value_conc1(jj,neighbors[ii][0],neighbors[ii][1],neighbors[ii][2])-value_conc1(jj,i,j,k) )*C_diffusivity(jj);

  del_conc*=(diff_dt/(dx*dx));

  total_conc(i,j,k)+=del_conc;

  Array<double,1> temp_conc(value_conc2(Range::all(),i,j,k));
  Partitioning( total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc );
  return;
}*/

void phase_field::Diffusion_solver(const int i,const int j,const int k,const real diff_dt) {
    //Array<int,2> neighbors(3,6,fortranArray);
    int neighbors[6][3];

    for(int ii=0; ii<6; ii++) {
        neighbors[ii][0]=i;
        neighbors[ii][1]=j;
        neighbors[ii][2]=k;
    }

    neighbors[0][0]=Get_neighbor(i,ix,-1);
    neighbors[1][0]=Get_neighbor(i,ix,1 );
    neighbors[2][1]=Get_neighbor(j,iy,-1);
    neighbors[3][1]=Get_neighbor(j,iy,1 );
    neighbors[4][2]=Get_neighbor(k,iz,-1);
    neighbors[5][2]=Get_neighbor(k,iz,1 );
    real del_conc=0.;

    for(int ii=0; ii<2*dimensions; ii++)
        for(int jj=1; jj<=max_num_phase; jj++)
            del_conc+=sqrt( value_phase(jj,neighbors[ii][0],neighbors[ii][1],neighbors[ii][2])*value_phase(jj,i,j,k) )*
                      ( value_conc1(jj,neighbors[ii][0],neighbors[ii][1],neighbors[ii][2])-value_conc1(jj,i,j,k) )*C_diffusivity(jj);

    del_conc*=(diff_dt/(dx*dx));
    total_conc(i,j,k)+=del_conc;
    Partitioning( total_conc(i,j,k),&value_phase(1,i,j,k),&value_conc2(1,i,j,k) );
    //Array<double,1> temp_conc(value_conc2(Range::all(),i,j,k));
    //Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);
    return;
}

real phase_field::Driving_force(const int grainii,const int phaseii,const int grainjj,const int phasejj,const real concii,const real concjj) const {
    real df=0.;

//     if(grainii<max_num_grains_per_phase)
//         df-=2e-12;//rex_energy(grainii);
//     else if(grainjj<max_num_grains_per_phase)
//         df+=2e-12;//rex_energy(grainjj);
    if(phaseii!=phasejj)
        //int i=max(phaseii,phasejj);
        //int j=min(phaseii,phasejj);
        df+=thermodynamic.para_entropy(phaseii,phasejj)*(concii-thermodynamic.para_equilibrium_conc(phaseii,phasejj))/thermodynamic.para_equilibrium_dc_dT(phaseii,phasejj);

    return df;
}

void phase_field::Partitioning(const real local_conc,const real* phase_array,real* conc_array) {
    for(int i=1; i<=max_num_phase; i++) {
        real t1=0.,t2=0.;
        const real *frac_ptr=phase_array;

        for (int j=1; j<=max_num_phase; j++) {
            t1+=*frac_ptr*thermodynamic.para_partition_k(i,j);
            t2+=*frac_ptr*thermodynamic.para_partition_c(i,j);
            frac_ptr++;
        }

        *conc_array=(local_conc-t2)/t1;
        conc_array++;
    }

    return;
}

// void phase_field::Partitioning(const double local_conc,const Array<double,1>& phase_array,Array<double,1> & conc_array)
// {
//     for(int i=1;i<=max_num_phase;i++)
//     {
//         double t1=0.,t2=0.;
//         for (int j=1;j<=max_num_phase;j++)
//         {
//             t1+=phase_array(j)*thermodynamic.para_partition_k(i,j);
//             t2+=phase_array(j)*thermodynamic.para_partition_c(i,j);
//             //cout<<t1<<" "<<thermodynamic.para_partition_k(i,j)<<" ";
//         }
//         conc_array(i)=(local_conc-t2)/t1;
// 
//     }
//     //cout<<endl;
//     //cout<<phase_array<<endl;
//     //cout<<conc_array<<endl;
//     return;
// }
void phase_field::Sort_xyz(Array<real,1>& value,Array<int,1>& grain,const int N) {
    for(int i=1; i<=N; i++) {
        TinyVector<int,1> max_index=maxIndex(value(Range(i,N)))+(i-1);
        std::swap(value(i),value(max_index));
        std::swap(grain(i),grain(max_index));
    }
}

