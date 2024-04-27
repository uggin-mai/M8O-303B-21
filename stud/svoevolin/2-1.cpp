#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include "matrix.h"

using namespace std;

using double_n = vector <double (*)(vector <double>)>;
using double_nn = vector <vector <double (*)(vector <double>)>>;
int iter = 0;

double dichotomy(function <double(double)> f, double l, double r, double eps)
{
    while (abs(l - r) > eps)
    {
        double nl = l + 0.3 * (r - l);
        double nr = r - 0.3 * (r - l);
        if (f(nl) < f(nr))
            l = nl;
        else
            r = nr;
    }
    return (l + r) / 2;
}

double newton(double (*f)(double), double (*df)(double), double x, double eps)
{
    double xk = x - f(x) / df(x);
    iter = 0;
    while (abs(xk - x) > eps)
    {
        x = xk;
        xk = x - f(x) / df(x);
        iter++;
    }
    return xk;
}

double iteration(double (*phi)(double), double (*dphi)(double), double l, double r, double eps)
{
    double x = (l + r) / 2;
    double xk = phi(x);
    iter = 0;
    double q = dphi(dichotomy(dphi, l, r, eps));
    while ((q * abs(xk - x) / (1 - q)) > eps)
    {
        x = xk;
        xk = phi(x);
        iter++;
    }
    return xk;
}

double f(double x)
{
    return tan(x) - 5. * pow(x, 2) + 1;
}

double df(double x)
{
    return (1. / pow(cos(x), 2)) - 10. * x;
}

double phi(double x)
{
    return atan(5. * pow(x, 2) - 1.);
}

double dphi(double x)
{
    return (10*x) / (pow((5 * pow(x, 2) - 1), 2) + 1);
}

int main()
{
    setlocale(LC_ALL, "Rus");
    ofstream fout("answer2-1.txt");
    fout.precision(4);
    fout << fixed;
    double ans = newton(f, df, 1.5, 0.00001);
    fout << "Ньютон: x=" << ans << ", f(x)=" << f(ans) << '\n';
    fout << "Количество итераций: " << iter << '\n';

    ans = iteration(phi, dphi, 1, 2, 0.00001);
    fout << "Простые итерации: x=" << ans << ", f(x)=" << f(ans) << '\n';
    fout << "Количество итераций: " << iter << '\n';
}