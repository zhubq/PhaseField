#ifndef _LINEAR_FUNCTOR_H
#define _LINEAR_FUNCTOR_H
#include "real.h"
class linear_functor {
public:
    linear_functor(real a=0.,real b=0.):a_(a),b_(b) {}
    void Reset(real a,real b) {
        a_=a;
        b_=b;
    }
    real operator ()(real T) {
        return (a_*T+b_);
    }
private:
    real a_,b_;
};
#endif