#include <math.h>
#include "newton.h"

using namespace std;



pair<double,double>  Newton::result(double eps0){
    double x01 = 1;
    double x02 = 1;
    double x1 = x01 - det_A1(x01,x02) / det_I(x01,x02);
    double x2 = x02 - det_A2(x01,x02) / det_I(x01,x02);
    while (eps(x1,x2,x01,x02) >  eps0){
        x01 = x1;
        x02 = x2;
        x1 = x01 - det_A1(x01,x02) / det_I(x01,x02);
        x2 = x02 - det_A2(x01,x02) / det_I(x01,x02);
    }

    pair<double,double>  res= {x1,x2};
    return res;



}
double  Newton::f1(double x1,double x2  ){
    return  3*pow(x1,2) - x1 + pow(x2,2) - 1;
}
double  Newton::f1_dx1(double x1,double x2  ){
    return  6*x1 - 1 ;
}
double Newton:: f1_dx2(double x1,double x2  ){
    return 2*x2;
}
double  Newton::f2(double x1,double x2  ){
    return x2 - tan(x1);
}
double  Newton::f2_dx1(double x1,double x2  ){
    return -1/pow(cos(x1),2);
}

double  Newton::f2_dx2(double x1,double x2  ){
    return 1;
}
double  Newton::det_I(double x1,double x2){
    double det = f1_dx1(x1,x2) * f2_dx2(x1,x2) - f2_dx1(x1,x2) * f1_dx2(x1,x2);
    return det;
}
double  Newton::det_A1(double x1,double x2){
    double det = f1(x1,x2) * f2_dx2(x1,x2) - f2(x1,x2) * f1_dx2(x1,x2);
    return det;
}
double  Newton::det_A2(double x1,double x2){
    double det = f1_dx1(x1,x2) * f2(x1,x2) - f2_dx1(x1,x2) * f1(x1,x2);
    return det;
}

double  Newton::eps(double x1,double x01,double x2,double x02){
//    return sqrt(pow(x2 - x02 ,2) + pow(x1 - x01 ,2));
    return max(x2 - x02 , x1 - x01);

}


