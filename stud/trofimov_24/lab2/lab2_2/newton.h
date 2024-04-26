

#ifndef LAB2_2_NEWTON_H
#define LAB2_2_NEWTON_H

#include <iostream>

using namespace std;
class Newton
{
public:
    pair<double,double>  result(double eps0);
private:
    double f1(double x1,double x2  );
    double f1_dx1(double x1,double x2  );
    double f1_dx2(double x1,double x2  );
    double f2(double x1,double x2  );
    double f2_dx1(double x1,double x2);
    double f2_dx2(double x1,double x2 );
    double det_I(double x1,double x2);
    double det_A1(double x1,double x2);
    double det_A2(double x1,double x2);
    double eps(double x1,double x01,double x2,double x02);
};

#endif //LAB2_2_NEWTON_H
