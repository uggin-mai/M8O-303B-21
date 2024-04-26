
#include "newton.h"

using namespace std;

double f(double x ){
    return pow(x,6 )- 5*x - 2;
}
double f_dev1(double x ){
    return 6*pow(x,5 )- 5;
}


double Newton(double x0,double eps){

    double x =  x0 - f(x0) / f_dev1(x0);
    double diff =  abs(x - x0);
    x0 = x;
    int k = 0;
    while (  diff >= eps ){
        x =  x0 - f(x0) / f_dev1(x0);
        diff =  abs(x - x0);
        x0 = x;
        k++;

    }
    return  x;

}
