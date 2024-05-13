#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f(double x) {
    return pow(3, x) - 5*x*x + 1;
}

double derivative(double x) {
    return log(3)*pow(3, x)-10*x;
}

double phi(double x){
    return sqrt((pow(3, x)+1)/5);
}

double simple_iterations(double x_0, double eps, int &iter){
    double x_i, x;
    x = x_0;
    for (int i = 0; i < 10000; i++) {
        x_i = phi(x);
        if (fabs(x_i - x) < eps)
            return x_i;
        x = x_i;
        iter++;
    }
    return x;
}

double newton(double x_0, double eps, int &iter){
    double x;
    x = x_0;
    while (fabs(f(x)) > eps && iter < 1000) {
        x = x - f(x) / derivative(x);
        iter++;
    }
    return x;
}

int main(){
    double x_0 = 0.5;
    double eps = 0.0001;
    int iterations = 0;
    ofstream fout("output.txt");
    double answer = simple_iterations(x_0, eps, iterations);
    fout << "Метод простых итераций:" << endl;
    fout << "x = " << answer << endl;
    fout << "Число итераций: " << iterations << endl;
    iterations = 0;
    answer = newton(x_0, eps, iterations);
    fout << "Метод Ньютона:" << endl;
    fout << "x = " << answer << endl;
    fout << "Число итераций: " << iterations << endl;
    return 0;
}