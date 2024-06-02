#include "newton.h"
using namespace std;
double f_iter(double x ){
    return 2/ (pow(x,5 )- 5);
}

double Iter(double x0,double eps){
    double x = f_iter(x0);
    double diff =  abs(x - x0);
    x0 = x;
    int k = 0;
    while (  diff >= eps ){
        x =  f_iter(x0);
        diff =  abs(x - x0);
        x0 = x;
        k++;

    }
    return  x;


}