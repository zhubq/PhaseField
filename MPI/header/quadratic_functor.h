#ifndef _QUADRATIC_FUNCTOR_H
#define _QUADRATIC_FUNCTOR_H

class quadratic_functor{
public:
    quadratic_functor(double a=0.,double b=0.,double c=0.):a_(a),b_(b),c_(c){}
    void Reset(double a,double b,double c){a_=a;b_=b;c_=c;}
    double operator ()(double x){
        return(a_*x*x+b_*x+c_);
    }
private:
    double a_,b_,c_;
};

#endif