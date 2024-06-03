#include <iostream>
#include <cmath>

using namespace std;

double module(double x) {
    return pow(pow(x,2),0.5);
}

double f1_x1_x2(double x1, double x2) {
    return pow(x1/2,2) + pow(x2/(2/2),2) - 1;
}

double f2_x1_x2(double x1, double x2) {
    return 2 * x2 - exp(x1) - x1;
}

double df1_dx1(double x1) {
    return x1/2;
}

double df1_dx2(double x2) {
    return 2*x2;
}

double df2_dx1(double x1) {
    return -exp(x1) - 1;
}

double df2_dx2() {
    return 2;
}

void iterations(double precision, double x1, double x2) {
    double x1_new;
    double x2_new;
    double x1_old = x1 + 1;
    double x2_old = x2 + 1;
    int counter = 0;
    while (module((x1_old - x1) > precision) || (module(x2_old - x2) > precision)) {
        if((pow(x2,2) > 1))
            break;
        x1_new = 2*pow((1-pow(x2,2)), 0.5);
        x2_new = (exp(x1) + x1)/2;
        x1_old = x1;
        x2_old = x2;
        x1 = x1_new;
        x2 = x2_new;
        counter++;
    }
    cout << "Iterations answer: " << "x1 = " << x1 << ", x2 = " << x2 << ", iterations = " << counter << endl;
    return;
}

void newton_method(double precision, double x1, double x2) {
    double dx1 = 1;
    double dx2 = 1;
    int counter = 0;
    while ((module(dx1) > precision || module(dx2) > precision)){
        double det = df1_dx1(x2) * df2_dx2() - df1_dx2(x1) * df2_dx1(x1);
        dx1 = (f1_x1_x2(x1, x2) * df2_dx2() - f2_x1_x2(x1, x2) * df1_dx2(x1)) / det;
        dx2 = (f2_x1_x2(x1, x2) * df1_dx1(x2) - f1_x1_x2(x1, x2) * df2_dx1(x1)) / det;

        x1 -= dx1;
        x2 -= dx2;
        counter++;
    }
    cout << "Newton answer:     " << "x1 = " << x1 << ", x2 = " << x2 << ", iterations = " << counter << endl;
    return;
}

int main() {
    double precision = 0.01;
    double x1 = 0.1, x2 = 0.1;
    iterations(precision, x1, x2);
    newton_method(precision, x1, x2);
    return 0;
}