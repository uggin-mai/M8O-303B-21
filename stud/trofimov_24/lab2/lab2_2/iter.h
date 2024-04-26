#ifndef LAB2_2_ITER_H
#define LAB2_2_ITER_H
#include <iostream>
using namespace std;
class Iter{
private:
    double  phi1(double x1,double x2  );
    double  phi2(double x1,double x2  );
    double eps(double x1,double x01,double x2,double x02);
public:
    pair<double,double> result(double eps0);

};
#endif //LAB2_2_ITER_H
