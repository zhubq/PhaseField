#ifndef _QUADRATIC_FUNCTOR_H
#define _QUADRATIC_FUNCTOR_H
#include "real.h"
class quadratic_functor {
public:
    quadratic_functor(real a=0.,real b=0.,real c=0.):a_(a),b_(b),c_(c) {}
    void Reset(real a,real b,real c) {
        a_=a;
        b_=b;
        c_=c;
    }
    real operator ()(real x) {
        return(a_*x*x+b_*x+c_);
    }
private:
    real a_,b_,c_;
};

#endif