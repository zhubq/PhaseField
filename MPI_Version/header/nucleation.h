#ifndef _NUCLEATION_H
#define _NUCLEATION_H
#include "phase_field.h"
#include "constant.h"
#include "thermo_data.h"
#include "random/uniform.h"
#include<algorithm>
class nucleation{
public:
    nucleation(){}
    void Nucleate1in2(phase_field &,const int &);
    
};


#endif