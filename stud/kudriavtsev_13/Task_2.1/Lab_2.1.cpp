#include <iostream>
#include <cmath>

using namespace std;

const double current_value = 3;

double f_x(double x) {
    return (log(x+1) - 2*x + 0.5);
}

double df_dx(double x) {
    return (1 / (x+1) - 2);
}

double g_x(double x){
    return (log(x+1)/2 + 0.25);
}

double module(double x) {
    return pow(pow(x,2),0.5);
}

double iterations(double x0, double precision) {
    double x = g_x(x0);
    double diff =  abs(x - x0);
    x0 = x;
    int k = 0;
    while (diff >= precision) {
        x =  g_x(x0);
        diff =  abs(x - x0);
        x0 = x;
        k++;
    }
    cout << endl;
    cout << "Iteration method answer: " << x << endl;
    cout << k << " iterations" << endl;
    return x;
}

double newton_method(double (*f_x)(double), double (*df_dx)(double), double precision) {
    double now, before = current_value;
    int iterations = 0;
    while (module((f_x(now))) > precision){
        now = before - f_x(before) / df_dx(before);
        iterations++;
        before = now;
    };
    cout << endl << "Newton's method answer: " << now << endl  << iterations << " iterations" << endl;
    return now;
}


int main() {
    double precision;
    cout << "Input precision: ";
    cin >> precision;
    iterations(current_value, precision);
    newton_method(f_x, df_dx, precision);
    return 0;
}