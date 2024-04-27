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

vector <double> newton_system(double_n f, double_nn df, vector <double> x, double eps)
{
    int n = f.size();
    matrix x0(n, 1);
    for (int i = 0; i < n; i++)
        x0[i][0] = x[i];
    matrix xk = x0;
    iter = 0;
    do
    {
        iter++;
        x0 = xk;
        for (int i = 0; i < n; i++)
            x[i] = x0[i][0];
        matrix a(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                a[i][j] = df[i][j](x);
        }
        matrix b(n, 1);
        for (int i = 0; i < n; i++)
            b[i][0] = -f[i](x);
        xk = x0 + solve_gauss(lu_decomposition(a), b);
    } while (abs(xk - x0) > eps);
    vector <double> res(n);
    for (int i = 0; i < n; i++)
        res[i] = xk[i][0];
    return res;
}

vector <double> iteration_system(double_n phi, double_nn dphi, vector <double> l, vector <double> r, double eps)
{
    double n = phi.size();
    vector <double> xk(n), x(n);
    for (int i = 0; i < n; i++)
        xk[i] = (l[i] + r[i]) / 2;
    vector <double> qn = xk;
    auto dphinx = [&](vector <double> x)
    {
        double res = 0;
        for (int j = 0; j < n; j++)
        {
            double tmp = 0;
            for (int k = 0; k < n; k++)
                tmp += abs(dphi[j][k](x));
            res = max(res, tmp);
        }
        return res;
    };
    for (int i = 0; i < n; i++)
    {
        auto dphix = [&](double x)
        {
            vector <double> xn = qn;
            xn[i] = x;
            return dphinx(xn);
        };
        qn[i] = dichotomy(dphix, l[i], r[i], eps);
    }
    double q = dphinx(qn);
    double epsk = 0;
    iter = 0;
    do
    {
        iter++;
        x = xk;
        for (int i = 0; i < n; i++)
            xk[i] = phi[i](x);
        epsk = 0;
        for (int i = 0; i < n; i++)
        {
            double s = abs(x[i] - xk[i]);
            epsk = max(epsk, s);
        }
    } while (q * epsk / (1 - q) > eps);
    return xk;
}

double f1(vector <double> x)
{
    return x[0] * x[0] - 2. * log10(x[1]) - 1.;
}

double f2(vector <double> x)
{
    return x[0] * x[0] - 2. * x[0] * x[1] + 2.;
}

double_n fn = { f1, f2 };

double df11(vector <double> x)
{
    return 2. * x[0];
}
double df12(vector <double> x)
{
    return - 2. / (x[1] * log(10));
}

double df21(vector <double> x)
{
    return 2 * x[0] - 2 * x[1];
}

double df22(vector <double> x)
{
    return - 2 * x[0];
}

double_nn dfn = { {df11, df12}, {df21, df22} };

double phi1(vector <double> x)
{
    return sqrt(2. * log10(x[1]) + 1);
}

double phi2(vector <double> x)
{
    return (pow(x[0], 2) + 2) / (2 * x[0]);
}

double_n phin = { phi1, phi2 };

double dphi11(vector <double> x)
{
    return 0;
}

double dphi12(vector <double> x)
{
    return pow((x[1] * sqrt(log(10))) * sqrt(2*log(x[1]) + log(10)), -1);
}

double dphi21(vector <double> x)
{
    return pow(2, -1) - pow(x[0], -2);
}

double dphi22(vector <double> x)
{
    return 0;
}

double_nn dphin = { {dphi11, dphi12}, {dphi21, dphi22} };

int main()
{
    setlocale(LC_ALL, "Rus");
    ofstream fout("answer2-2.txt");
    fout.precision(4);
    fout << fixed;
    vector <double> x = newton_system(fn, dfn, { 1, 1 }, 0.00001);
    fout << "Система методом Ньютона: x1=" << x[0] << ", x2=" << x[1] << ", f1(x)=" << f1(x) << ", f2(x)=" << f2(x) << '\n';
    fout << "Количество итераций: " << iter << '\n';

    x = iteration_system(phin, dphin, { 1, 1 }, { 2, 2 }, 0.00001);
    fout << "Система методом простых итераций: x1=" << x[0] << ", x2=" << x[1] << ", f1(x)=" << f1(x) << ", f2(x)=" << f2(x) << '\n';
    fout << "Количество итераций: " << iter << '\n';
}
