#ifndef _NUCLEATION_H
#define _NUCLEATION_H
#include "phase_field.h"
#include "constant.h"
#include "thermo_data.h"
#include "random/uniform.h"
#include<algorithm>
#include "real.h"
/*
 * A class to manipulate nucleation 
 */
class nucleation {
public:
    nucleation() {}
    void Nucleate1in2(phase_field &);//Nucleation of Phase 1 in Phase 2;

};


#endif