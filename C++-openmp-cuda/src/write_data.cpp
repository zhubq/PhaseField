#include "phase_field.h"
#include <fstream>
void phase_field::Write_data(bool forced ) {
    static real last_temperature=0.;
    static int micro_count=0;

    if ( forced||(abs(time-timestep_output*floor(time/timestep_output+0.5))<0.5*pfm_dt)||
            abs(temperature-last_temperature)>10.) {
        kinetics<<setw(8)<<setprecision(4)<<time<<" "<<temperature<<" ";
        Array<float,1> vol_fraction(max_num_phase,fortranArray);
        real domain_vol=ix*iy*iz;

        for(int i=1; i<=max_num_phase; i++) {
            vol_fraction(i)=sum( value_phase(i,Range::all(),Range::all(),Range::all()) )/domain_vol;
            kinetics<<setprecision(4)<<vol_fraction(i)<<" ";
        }

        kinetics<<sum( where(index_grain1<max_num_grains_per_phase, value_grain1, 0.) )/domain_vol<<endl;
        cout<<"time:"<<time<<" "<<"dt: "<<pfm_dt<<" "<<diff_dt<<endl;

        if (Is_micro_out) {
            ostringstream filename;
            filename<<"field";
            filename.width(7);
            filename.fill('0');
            filename<<int(time*100)<<"_"<<micro_count<<".vti";
            ofstream microstructure( (directory+'/'+filename.str()).c_str());
            // void vti_write(ofstream &,Array<int,4>&,Array<real,3>&,real);
//             cout<<index_grain1.lbound()<<endl;
//             cout<<index_grain1.ubound()<<endl;
            Vti_write(microstructure,index_grain1,total_conc,dx);
            microstructure.close();
            micro_count++;
        }

        last_temperature=temperature;
    }

    return;
}

void phase_field::Vti_write(ofstream & file,Array<int,4>& index,Array<real,3>& conc,real dx) {
    file<<"<?xml version=\"1.0\"?>"<<'\n'
        <<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<"\n"
        <<"<ImageData WholeExtent=\""<<index.lbound(1)<<' '<<index.ubound(1)<<' '
        <<index.lbound(2)<<' '<<index.ubound(2)<<' '
        <<index.lbound(3)<<' '<<index.ubound(3)<<'"'
        <<" Origin="<<"\"0 0 0\" Spacing=\""<<dx<<' '<<dx<<' '<<dx<<"\">"<<endl;
    file<<"<Piece Extent=\""<<index.lbound(1)<<' '<<index.ubound(1)<<' '
        <<index.lbound(2)<<' '<<index.ubound(2)<<' '
        <<index.lbound(3)<<' '<<index.ubound(3)<<"\">"<<endl;
    file<<"<PointData Scalars=\"energy\">"<<endl;
    file<<"<DataArray type=\"Int32\" Name=\"korn\" format=\"ascii\">"<<endl;

    for(int k=index.lbound(3); k<=index.ubound(3); k++)
        for(int j=index.lbound(2); j<=index.ubound(2); j++)
            for(int i=index.lbound(1); i<=index.ubound(1); i++)

//                  cout<<i<<" "<<j<<" "<<k<<endl;
                if(index(2,i,j,k)>0)
                    file<<0<<" ";

                else
                    file<<Phase_of_grain(index(1,i,j,k))<<" ";

    // file<<(index(1,i,j,k))<<" ";
    file<<"\n</DataArray>"<<endl;
    file<<"<DataArray type=\"Float32\" Name=\"conc\" format=\"ascii\">"<<endl;

    for(int k=conc.lbound(2); k<=conc.ubound(2); k++)
        for(int j=conc.lbound(1); j<=conc.ubound(1); j++)
            for(int i=conc.lbound(0); i<=conc.ubound(0); i++) {
//                  cout<<conc.lbound()<<endl<<conc.ubound()<<endl;
                file<<setprecision(3)<<total_conc(i,j,k)<<" ";
            }

    file<<"</DataArray>"<<endl
        <<"</PointData>"<<endl
        <<"</Piece>"<<endl
        <<"</ImageData>"<<endl
        <<"</VTKFile>"<<endl;
}
