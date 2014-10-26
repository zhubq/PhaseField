#include "phase_field.h"
#include "random/uniform.h"

void phase_field::Construct(bool c,char a[]) {
    if (c) {
        cout<<"using Voronoi tessellation to construct the structure."<<endl;
        Array<int, 3> grain_index(ix,iy,iz,fortranArray);
        Tessellation(grain_index,num_grain[1],reduction);

        for(int k=1; k<=iz; k++)
            for(int j=1; j<=iy; j++)
                for(int i=1; i<=ix; i++)
                    index_grain1(1,i,j,k)=grain_index(i,j,k)+200000;

    } else {
        cout<<"read initial microstructure from file "<<a<<endl;
        ifstream file(a);

        for(int k=1; k<=iz; k++)
            for(int j=1; j<=iy; j++)
                for(int i=1; i<=ix; i++)
                    file>>index_grain1(1,i,j,k);

        file.close();
    }

    value_grain1(1,Range::all(),Range::all(),Range::all())=1.0;
    #pragma omp parallel for collapse(2)

    for(int k=1; k<=iz; k++)
        for(int j=1; j<=iy; j++)
            for (int i=1; i<=ix; i++) {
                value_phase(Phase_of_grain(index_grain1(1,i,j,k)),i,j,k)=1.0;
            }

    for(int i=0; i<max_num_phase; i++)
        value_conc1(i+1,Range::all(),Range::all(),Range::all())=init_conc[i];

    cout<<"start initial iteration to generate diffuse interfaces"<<endl;

    for(int n=1; n<=num_init_interation; n++) {
        #pragma omp parallel for collapse(2)
        for(int k=1; k<=iz; k++)
            for(int j=1; j<=iy; j++)
                for(int i=1; i<=ix; i++) {
                    Iso_growth(i,j,k);
                }

        cycleArrays(index_grain1,index_grain2);
        cycleArrays(value_grain1,value_grain2);
    }

    Update_parameter();
    #pragma omp for collapse(2)
    for(int k=1; k<=iz; k++)
        for(int j=1; j<=iy; j++)
            for(int i=1; i<=ix; i++) {
                total_conc(i,j,k)=dot( value_phase(Range::all(),i,j,k),value_conc1(Range::all(),i,j,k) );
                Partitioning(total_conc(i,j,k),&value_phase(1,i,j,k),&value_conc1(1,i,j,k));
	        //Array<double,1> temp_conc(value_conc1(Range::all(),i,j,k));
                //Partitioning(total_conc(i,j,k),value_phase(Range::all(),i,j,k),temp_conc);
            }

    if(num_init_interation>0)cout<<"******Diffuse Interface Generated******"<<endl;//cout<<pi*pi/(2.*eta*eta)*(dx*dx)<<endl;
}


void   phase_field::Tessellation(Array<int,3> & grain_index,int num, real reduction) {
    ranlib::UniformClosed<real> uniform_numbers;
    uniform_numbers.seed(num);
    Array<real,2> point_xyz(3,num,fortranArray);

    for(Array<real,2>::iterator i=point_xyz.begin(); i!=point_xyz.end(); i++) {
        *i=uniform_numbers.random();
    }

    //cout<<point_xyz<<endl;
    //cout<<uniform_numbers.random()<<endl;
    point_xyz(1,Range::all())*=(grain_index.ubound(0)-grain_index.lbound(0));
    point_xyz(1,Range::all())+=grain_index.lbound(0);
    point_xyz(2,Range::all())*=(grain_index.ubound(1)-grain_index.lbound(1));
    point_xyz(2,Range::all())+=grain_index.lbound(1);
    point_xyz(3,Range::all())*=(grain_index.ubound(2)-grain_index.lbound(2));
    point_xyz(3,Range::all())+=grain_index.lbound(2);

    TinyVector<int,3> lbound=grain_index.lbound(), ubound=grain_index.ubound(),shape=grain_index.shape();
    
    #pragma omp parallel for shared(grain_index,num)
    for(int k=lbound(2); k<=ubound(2); k++)
        for(int j=lbound(1); j<=ubound(1); j++)
            for(int i=lbound(0); i<=ubound(0); i++) {
                int index=0;
                real temp1=dot(shape/(1.-reduction),shape/(1.-reduction));
                real temp2;

                for(int n=1; n<=num; n++) {
                    temp2=pow2( ( abs(point_xyz(1,n)-i)<0.5*shape(0)?abs(point_xyz(1,n)-i):shape(0)-abs(point_xyz(1,n)-i) )/(1.-reduction) )
                          +pow2( ( abs(point_xyz(2,n)-j)<0.5*shape(1)?abs(point_xyz(2,n)-j):shape(1)-abs(point_xyz(2,n)-j) )*(1.-reduction) )
                          +pow2( abs(point_xyz(3,n)-k)<0.5*shape(2)?abs(point_xyz(3,n)-k):shape(2)-abs(point_xyz(3,n)-k) );

                    if(temp1>temp2) {
                        temp1=temp2;
                        index=n;
                    }
                }

                grain_index(i,j,k)=index;
            }
}

static void Select_bands(Array<int,3> & grain_index,int num,real fraction,int spacing) {
    int num_pixel=grain_index.size()*fraction;
    Array<int,1> grain_volume(num,fortranArray);

    for(int i=1; i<=num; i++)
        grain_volume(i)=count(grain_index==i);

    Array<int,1> Is_InBands(num,fortranArray);
    
    #pragma omp for
    for (int k=grain_index.lbound(2); k<=grain_index.ubound(2); k++)
        for(int j=grain_index.lbound(1)+1; j<grain_index.ubound(1); j+=spacing)
            for(int i=grain_index(0); i<=grain_index.ubound(0); i++) {
                ++Is_InBands(grain_index(i-1,j,k));
                ++Is_InBands(grain_index(i,j,k));
                ++Is_InBands(grain_index(i+1,j,k));
            }
}

void phase_field::Iso_growth(const int i,const int j,const int k) {
    
    const static int max_num_nonzero_stencil=10;
    Array<int,1> index_grain_nonzero_stencil(max_num_nonzero_stencil,fortranArray);//assume there are at most 10 non-zero phase-fields on all stencil-points,7-stencil pattern is used in 3D.
    index_grain_nonzero_stencil=0;
    Array<int,1> index_phase_nonzero_stencil(max_num_nonzero_stencil,fortranArray);
    Array<real,2> value_grain_nonzero_stencil(max_num_nonzero_stencil,2,fortranArray);
    Array<real,1> del_value(max_num_nonzero_stencil,fortranArray);
    //Array<real,1> local_conc(value_conc1(Range::all(),i,j,k));
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

    for(int ii=1; ii<=max_num_coexist; ii++) {
        if(value_grain1(ii,i,j,k)>cutoff) {
            num_nonzero_stencil++;
            index_grain_nonzero_stencil(ii)=index_grain1(ii,i,j,k);
            value_grain_nonzero_stencil(ii,1)=value_grain1(ii,i,j,k);
        }
    }

    //num_nonzero_stencil=num_nonzero_local;
    for(int ii=1; ii<=2*dimensions; ii++) {
        for(int jj=1; jj<=max_num_coexist; jj++) {
            int temp_index=index_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));

            if(temp_index==0)break;

            real temp_value=value_grain1(jj,neighbors(1,ii),neighbors(2,ii),neighbors(3,ii));
            int kk;

            if(temp_value>cutoff) {
                for(kk=1; kk<=max_num_nonzero_stencil; kk++) {
                    if(temp_index==index_grain_nonzero_stencil(kk)) {
                        value_grain_nonzero_stencil(kk,2)+=temp_value;
                        break;
                    }
                }
            }

            if(kk>max_num_nonzero_stencil && num_nonzero_stencil<max_num_nonzero_stencil) {
                num_nonzero_stencil+=1;
                index_grain_nonzero_stencil(num_nonzero_stencil)=temp_index;
                value_grain_nonzero_stencil(num_nonzero_stencil,2)+=temp_value;
            }
        }
    }

    if(num_nonzero_stencil==1) {
      for(int ii=1; ii<=num_nonzero_stencil;ii++) {
        index_grain2(ii,i,j,k)=index_grain1(ii,i,j,k);
        value_grain2(ii,i,j,k)=value_grain1(ii,i,j,k);
      }
        return;
    }

//---------------------------------------------------------------------------------------------
//=============================================================================================
    /*======now,three arrays to be used for calculation:
     * index_grain_nonzero_stencil,index_phase_nonzero_stencil,value_grain_nonzero_stencil.
     *
    */
    del_value=0.;
    real factor1=1./(dx*dx),
         factor2=pi*pi/(2.*eta*eta);
    real del_value_ij;

    for(int ii=1; ii<=num_nonzero_stencil; ii++) {
        for(int jj=1; (jj<=num_nonzero_stencil); jj++) {
            if (ii==jj)continue;

            del_value_ij=( value_grain_nonzero_stencil(ii,2)*value_grain_nonzero_stencil(jj,1)-
                           value_grain_nonzero_stencil(jj,2)*value_grain_nonzero_stencil(ii,1) )*factor1
                         +(value_grain_nonzero_stencil(ii,1)-value_grain_nonzero_stencil(jj,1))*factor2;
            del_value(ii)+=del_value_ij;
        }

        del_value(ii)*=0.5/dimensions*dx*dx;
        del_value(ii)+=value_grain_nonzero_stencil(ii,1);

        if(del_value(ii)>1.0)del_value(ii)=1.0;

        else if(del_value(ii)<cutoff)del_value(ii)=0.;
    }      
    
    Sort_xyz(del_value,index_grain_nonzero_stencil,max_num_nonzero_stencil);
    real sum_del_value=0.;
    value_phase(Range::all(),i,j,k)=0.;
    
    for(int ii=1; ii<=max_num_coexist; ii++) {
        if (del_value(ii)>cutoff) {
            index_grain2(ii,i,j,k)=index_grain_nonzero_stencil(ii);
            value_grain2(ii,i,j,k)=del_value(ii);
            int index=Phase_of_grain( index_grain_nonzero_stencil(ii) );
            value_phase(index,i,j,k)+=del_value(ii);
            sum_del_value+=del_value(ii);
        }
        else {
	  index_grain2(ii,i,j,k)=0;
	  value_grain2(ii,i,j,k)=0.;
	}
    }
    
    for (int ii=1;ii<=max_num_coexist;ii++)
        value_grain2(ii,i,j,k)*=1./sum_del_value;
   
    for(int ii=1;ii<=max_num_phase;ii++)
        value_phase(ii,i,j,k)*=1./sum_del_value;
    return;
}
