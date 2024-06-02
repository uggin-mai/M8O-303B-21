#include <math.h>
#include "iter.h"

using namespace std;




double Iter::phi1(double x1, double x2) {
    return  ((1 - pow(x2,2))/(3*x1 - 1) + x1) / 2;
}

double Iter::phi2(double x1, double x2) {
    return tan(x1);
}

double Iter::eps(double x1,double x01,double x2,double x02){
    return max(abs(x2 - x02) , abs(x1 - x01));

}
pair<double,double> Iter::result(double eps0){
    double x01 = 10 ;
    double x02 =10;

    double x1 = phi1(x01,x02);
    double x2 = phi2(x01,x02);

    while ( eps(x1,x01,x2,x02)> eps0){

        x01 = x1;
        x02 = x2;
        x1  = phi1(x01,x02);
        x2  = phi2(x01,x02);

    }

    pair<double,double>  res= {x1,x2};
    return res;



}

