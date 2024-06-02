#include <iostream>
#include "iter.h"
#include "newton.h"




int main() {

    double x0;
    cout << "init x0:";
    cin >> x0;

    double precision;
    cout << "init  precision:";
    cin >> precision;

    double  res1 = Newton(x0,precision);
    std::cout <<"Newton: " <<res1 << std::endl;

    std::cout <<"f(x) =  " <<f(res1) << std::endl;


    double  res2 = Iter(x0,precision);
    std::cout <<"Iter: " <<res2 << std::endl;
    std::cout <<"f(x) =  " <<f(res1) << std::endl;

    return 0;
}
