#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


double F(double x) {
    return  sin(x) - pow(x,2) + 1; //sin x - x^2 + 1
}

double Diff_F(double x) {
    return cos(x) - 2*x; //cos x - 2x
}

double Diff2_F(double x) {
    return -sin(x) - 2; //-sin x - 2
}

double Phi(double x) {
    return pow(1+sin(x), 0.5); //sqrt(1+sin x)
}

double Diff_Phi(double x) {
    return cos(x)/pow(1+sin(x), 0.5);//cos(x)/sqrt(1+sin x)
}

double Iterations_method(double x0, double eps, int &i) {
    int k = 0;
    double q = abs(Diff_Phi(0.75)), x1 = x0-eps-1;
    do{
        x0 = x1;
        x1 = Phi(x1);
        i += 1;
    }while (q / (1 - q) * abs(x1 - x0) > eps);
    return x1;
}
double Newton_method(double x0, double eps, int &i) {
    int k = 0;
    double x1 = x0-eps-1;;
    do{
        x0 = x1;
        x1 = x0 - F(x1) / Diff_F(x1);
        i += 1;
    }while (abs(x1 - x0) >= eps);
    return x1;
}

int main() {
    //Начальные приближения подобраны графически исходя из условий выполнения методов
    cout.precision(9);
    ofstream fout("answer.txt");
    ifstream fin("input.txt");
    double eps;
    fin >> eps;
    double X_iterations, x0_iterations = 1.375;
    int iterator_iterations = 0;
    X_iterations = Iterations_method(x0_iterations, eps, iterator_iterations);
    fout << "===Simple iterations method===\nIterations number: " << iterator_iterations << "\nRoot: " << to_string(X_iterations) << '\n';
    double X_Newton, x0_Newton = 1.7;
    int iterator_Newton = 0;
    X_Newton = Newton_method(x0_Newton, eps, iterator_Newton);
    fout << "===Newton method===\nIterations number: " << iterator_iterations << "\nRoot: " << to_string(X_iterations) << '\n';
}