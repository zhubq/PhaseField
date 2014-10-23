#ifndef _LINEAR_FUNCTOR_H
#define _LINEAR_FUNCTOR_H

class linear_functor{
public:
    linear_functor(double a=0.,double b=0.):a_(a),b_(b){}
    void Reset(double a,double b){a_=a;b_=b;}
    double operator ()(double T){
        return (a_*T+b_);
    }
private:
    double a_,b_;
};
#endif