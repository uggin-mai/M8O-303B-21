#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// Функция, корень которой мы ищем: ln(x+2) - x^2
double F(double x) {
    return log(x + 2) - pow(x, 2);
}

// Первая производная функции F
double Diff_F(double x) {
    return 1 / (x + 2) - 2 * x;
}

// Функция трансформации Phi для метода простых итераций
double Phi(double x) {
    return x + 0.1 * (log(x + 2) - pow(x, 2));
}

// Метод простых итераций
double Iterations_method(double x0, double eps, int& i) {
    double x1 = x0;
    do {
        x0 = x1;
        x1 = Phi(x1);
        i += 1;
    } while (abs(x1 - x0) > eps);
    return x1;
}

// Метод Ньютона
double Newton_method(double x0, double eps, int& i) {
    double x1 = x0;
    do {
        x0 = x1;
        x1 = x0 - F(x1) / Diff_F(x1);
        i += 1;
    } while (abs(x1 - x0) >= eps);
    return x1;
}

int main() {
    cout.precision(9);
    ofstream fout("answer.txt");
    double eps=0.000001;
    double X_iterations, x0_iterations = 1.0;  // начальное приближение
    int iterator_iterations = 0;
    X_iterations = Iterations_method(x0_iterations, eps, iterator_iterations);
    fout << "===Simple iterations method===\nIterations number: " << iterator_iterations << "\nRoot: " << to_string(X_iterations) << '\n';
    double X_Newton, x0_Newton = 1.0;  // начальное приближение
    int iterator_Newton = 0;
    X_Newton = Newton_method(x0_Newton, eps, iterator_Newton);
    fout << "===Newton method===\nIterations number: " << iterator_Newton << "\nRoot: " << to_string(X_Newton) << '\n';
}
