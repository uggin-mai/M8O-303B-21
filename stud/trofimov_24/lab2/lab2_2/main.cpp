#include <iostream>
#include "newton.h"
#include "iter.h"
#include "math.h"




double f(double x ){
    return pow(x,2) - 0.25 * x + 0.015625;
}

int main() {
    Newton nt;
    std::pair<double,double>   res1 = nt.result(0.01);
    std::cout <<"Newton: " <<res1.first<<" " <<res1.second<< std::endl;

    std::cout << std::endl;

    Iter it;
    std::pair<double,double>   res2 = it.result(0.01);
    std::cout <<"Iter: " <<res2.first<<" " <<res2.second<< std::endl;

//
//    double  res2 = Iter(0,0.01);
//    std::cout <<"Iter: " <<res2 << std::endl;
//    std::cout << f(0.084) << std::endl;

    return 0;
}
