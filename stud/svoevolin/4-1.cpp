#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>
#include "matrix.h"

using namespace std;
using equation = function <double(vector <double>)>;
using func = function <double(double)>;

double runge_romberg(double i1, double i2, double h1, double h2, double p)
{
    double k = h2 / h1;
    return (i1 - i2) / (pow(k, p) - 1);
}

vector <equation> make_system(equation dy, int k)
{
    vector <equation> res(k);
    for (int i = 0; i < k - 1; i++)
    {
        auto f = [=](vector <double> x) mutable
        {
            return x[i + 2];
        };
        res[i] = f;
    }
    res[k - 1] = dy;
    return res;
}

vector <vector <double>> explicit_euler(vector <equation> dy, vector <double> yk, double a, double b, double h)
{
    int n = dy.size();
    vector <vector <double>> res(n);
    for (int i = 0; i < n; i++)
        res[i].push_back(yk[i]);
    for (double x = a; x <= b - h; x += h)
    {
        vector <double> args;
        args.push_back(x);
        for (int i = 0; i < n; i++)
            args.push_back(res[i].back());
        for (int i = 0; i < n; i++)
            res[i].push_back(res[i].back() + h * dy[i](args));
    }
    return res;
}

vector <vector <double>> improved_euler(vector <equation> dy, vector <double> yk, double a, double b, double h)
{
    int n = dy.size();
    vector <vector <double>> tmp(n);
    vector <vector <double>> res(n);
    for (int i = 0; i < n; i++)
    {
        tmp[i].push_back(yk[i]);
        res[i].push_back(yk[i]);
    }
    int k = 0;
    for (double x = a; x <= b - h / 2; x += h / 2)
    {
        vector <double> args;
        args.push_back(x);
        for (int i = 0; i < n; i++)
            args.push_back(tmp[i].back());
        for (int i = 0; i < n; i++)
        {
            tmp[i].push_back(tmp[i][tmp[i].size() - 1 - k % 2] + 0.5 * (k % 2 + 1) * h * dy[i](args));
            if (k % 2)
                res[i].push_back(tmp[i].back());
        }
        k++;
    }
    return res;
}

vector <double> make_args(double x, const vector <double>& y, double add_x, double add_y)
{
    vector <double> res;
    res.push_back(x + add_x);
    for (int i = 0; i < y.size(); i++)
        res.push_back(y[i] + add_y);
    return res;
}

vector <vector <double>> runge_kutta(vector <equation> dy, vector <double> yk, double l, double r, double h)
{
    // для 4 порядка точности
    int p = 4;
    vector <double> a = { 0, 0, 0.5, 0.5, 1 };
    vector <vector <double>> b = { {}, {0}, {0, 0.5}, {0, 0, 0.5}, {0, 0, 0, 0.5} };
    vector <double> c = { 0, 1. / 6, 1. / 3, 1. / 3, 1. / 6 };
    //
    int n = dy.size();
    vector <vector <double>> res(n);
    for (int i = 0; i < n; i++)
        res[i].push_back(yk[i]);
    vector <double> K(p + 1);
    for (double x = l; x <= r - h; x += h)
    {
        vector <double> y;
        for (int idx = 0; idx < n; idx++)
            y.push_back(res[idx].back());
        for (int idx = 0; idx < n; idx++)
        {
            K[1] = h * dy[idx](make_args(x, y, 0, 0));
            for (int i = 2; i <= p; i++)
            {
                double add = 0;
                for (int j = 1; j <= i - 1; j++)
                    add += b[i][j] * K[j];
                K[i] = h * dy[idx](make_args(x, y, a[i] * h, add));
            }
            double delta = 0;
            for (int i = 1; i <= p; i++)
                delta += c[i] * K[i];
            res[idx].push_back(res[idx].back() + delta);
        }
    }
    return res;
}

vector <vector <double>> adams(vector <equation> dy, vector <double> yk, double a, double b, double h)
{
    int n = dy.size();
    vector <vector <double>> y = runge_kutta(dy, yk, a, b, h);
    int m = y[0].size();
    vector <vector <double>> res(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 4; j++)
            res[i].push_back(y[i][j]);
    }
    for (int k = 4; k < m; k++)
    {
        vector <vector <double>> args(4);
        for (int j = 0; j < 4; j++)
        {
            args[j].push_back(a + h * (k - j - 1));
            for (int i = 0; i < n; i++)
                args[j].push_back(res[i][k - j - 1]);
        }
        for (int i = 0; i < n; i++)
        {
            double delta = 55 * dy[i](args[0]) - 59 * dy[i](args[1]) + 37 * dy[i](args[2]) - 9 * dy[i](args[3]);
            res[i].push_back(res[i].back() + (h / 24) * delta);
        }
    }
    return res;
}

int main()
{
    setlocale(LC_ALL, "Rus");
    ofstream fout("answer4-1.txt");
    fout.precision(8);
    fout << fixed;

    auto ddy = [](vector <double> x)
        {
            return ( -1./2 * x[2] + 3./4 * x[1] ) / ( x[0] * (x[0] - 1) );
        };
    double h1 = 0.1, h2 = 0.05;
    vector <vector <double>> y1 = explicit_euler(make_system(ddy, 2), { 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h1);
    vector <vector <double>> y2 = improved_euler(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h1);
    vector <vector <double>> y3 = runge_kutta(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h1);
    vector <vector <double>> y4 = adams(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h1);
    fout << "Явный Эйлер" << "\t" << "Улучшенный" << "\t" << "Рунге-Кутта" << "\t" << "Адамс      " << "\t" << "Точный ответ" << endl;
    for (int i = 0; i < y1[0].size(); i++)
    {
        double x = 2 + h1 * i;
        fout << y1[0][i] << '\t' << y2[0][i] << '\t' << y3[0][i] << '\t' << y4[0][i] << '\t' << pow(abs(x), 3./2) << endl;
    }
    vector <vector <double>> y12 = explicit_euler(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h2);
    vector <vector <double>> y22 = improved_euler(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h2);
    vector <vector <double>> y32 = runge_kutta(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h2);
    vector <vector <double>> y42 = adams(make_system(ddy, 2),{ 2 * sqrt(2), 3. / 2 * sqrt(2) }, 2, 3, h2);
    fout << "\nПогрешности методом Рунге-Ромберга\n";
    fout << "Явный Эйлер" << "\t" << "Улучшенный" << "\t" << "Рунге-Кутта" << "\t" << "Адамс      " << endl;
    for (int i = 0; i < y1[0].size(); i++)
    {
        fout << runge_romberg(y1[0][i], y12[0][2 * i], h1, h2, 2) << '\t'
            << runge_romberg(y2[0][i], y22[0][2 * i], h1, h2, 4) << '\t'
            << runge_romberg(y3[0][i], y32[0][2 * i], h1, h2, 4) << '\t'
            << runge_romberg(y4[0][i], y42[0][2 * i], h1, h2, 4) << endl;
    }
}
