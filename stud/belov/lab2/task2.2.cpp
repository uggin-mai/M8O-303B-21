#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

// Функции системы уравнений
double f1(double x1, double x2) {
    return (x1 * x1 + 9) * x2 - 27;
}

double f2(double x1, double x2) {
    return (x1 - 1.5) * (x1 - 1.5) + (x2 - 1.5) * (x2 - 1.5) - 9;
}

// Якобиан системы
vector<vector<double>> jacobian(double x1, double x2) {
    vector<vector<double>> J(2, vector<double>(2));
    J[0][0] = 2 * x1 * x2;
    J[0][1] = x1 * x1 + 9;
    J[1][0] = 2 * (x1 - 1.5);
    J[1][1] = 2 * (x2 - 1.5);
    return J;
}

// Метод Ньютона
vector<double> newtonMethod(double x1, double x2, double tol) {
    vector<double> x = { x1, x2 };
    int iteration = 0;
    while (true) {
        vector<vector<double>> J = jacobian(x[0], x[1]);
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-6) break;

        vector<vector<double>> invJ(2, vector<double>(2));
        invJ[0][0] = J[1][1] / det;
        invJ[0][1] = -J[0][1] / det;
        invJ[1][0] = -J[1][0] / det;
        invJ[1][1] = J[0][0] / det;

        vector<double> f = { f1(x[0], x[1]), f2(x[0], x[1]) };
        vector<double> dx = { invJ[0][0] * f[0] + invJ[0][1] * f[1], invJ[1][0] * f[0] + invJ[1][1] * f[1] };
        x[0] -= dx[0];
        x[1] -= dx[1];
        iteration++;

        if (sqrt(dx[0] * dx[0] + dx[1] * dx[1]) < tol) break;
    }
    return x;
}

// Метод простой итерации
vector<double> simpleIteration(double x1, double tol) {
    double x2 = 27 / (x1 * x1 + 9);
    int iteration = 0;
    double x1_new;
    do {
        x1_new = x1;
        x1 = sqrt(9 - (x2 - 1.5) * (x2 - 1.5)) + 1.5;
        x2 = 27 / (x1 * x1 + 9);
        iteration++;
    } while (fabs(x1 - x1_new) > tol && iteration < 10000);
    return { x1, x2 };
}

int main() {
    ofstream fout("answer.txt");
    double x1_initial = 2.0;
    double x2_initial = 2.0;
    double tol = 1e-6;

    vector<double> result_newton = newtonMethod(x1_initial, x2_initial, tol);
    vector<double> result_si = simpleIteration(x1_initial, tol);

    fout << "Newton Method Result:\n";
    fout << "x1 = " << result_newton[0] << ", x2 = " << result_newton[1] << endl;
    fout << "Simple Iteration Result:\n";
    fout << "x1 = " << result_si[0] << ", x2 = " << result_si[1] << endl;
    fout.close();

    return 0;
}
