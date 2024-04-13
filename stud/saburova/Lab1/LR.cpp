#include <iostream>
#include <vector>
#include <ccomplex>
#include <fstream>
#include "matrix.h"

using namespace std;


pair<matrix, matrix> lu(matrix u, matrix& roots, bool root_flag)
{
    int n = u.n;
    matrix l(n, n);
    for (int k = 0; k < n; k++)
    {
        matrix prev = u;

        int swap_index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(prev[swap_index][k]) < abs(prev[i][k]))
                swap_index = i;
        }

        swap(u[k], u[swap_index]);
        swap(l[k], l[swap_index]);
        swap(prev[k], prev[swap_index]);
        if (root_flag)
            swap(roots[k], roots[swap_index]);

        for (int i = k + 1; i < n; i++)
        {
            double h = prev[i][k] / prev[k][k];
            l[i][k] = h;
            for (int j = k; j < n; j++)
                u[i][j] = prev[i][j] - h*prev[k][j];
        }
    }
    for (int i = 0; i < n; i++)
        l[i][i] = 1;
    return make_pair(l, u);
}

matrix direct_passage(matrix x_matr, matrix y_vect, bool flag)
{
    int n = x_matr.n;
    matrix res(n, 1);
    int d, first;
    if (flag) 
    {
        first = n-1;
        d = -1;
    } 
    else 
    { 
        first = 0;
        d = 1;
    }
    for (int i = first; i < n && i >= 0; i += d)
    {
        res[i][0] = y_vect[i][0];
        for (int j = 0; j < n; j++)
        {
            if (i != j)
                res[i][0] -= x_matr[i][j] * res[j][0];
        }
        res[i][0] = res[i][0] / x_matr[i][i];
    }
    return res;
}

matrix solve(pair <matrix, matrix> lu, matrix y_vect)
{
    matrix z = direct_passage(lu.first, y_vect, false);
    matrix x = direct_passage(lu.second, z, true);
    return x;
}

matrix inverse(matrix u, matrix roots)
{
    int n = u.n;
    matrix y_vect(n, 1);
    pair <matrix, matrix> lu_pair = lu(u, roots, true);
    matrix res(n, n);
    for (int i = 0; i < n; i++)
    {
        y_vect[max(i - 1, 0)][0] = 0;
        y_vect[i][0] = 1;
        matrix col = solve(lu_pair, y_vect);
        for (int j = 0; j < n; j++)
            res[j][i] = col[j][0];
    }
    return res;
}

double det(matrix u, matrix roots) 
{
    int n = u.n;
    pair <matrix, matrix> lu_pair = lu(u, roots, false);
    double res = 1;
    for (int i = 0; i < n; i++)
        res *= lu_pair.second[i][i];
    return res;
}


int main()
{
    cout << "If you want to change input data, open \"res\" directory and and change the data in the file of the corresponding task" << endl;
    cout << "Enter the code of task you want to check:" << endl;
    cout << "1) LU decomposition" << endl;
    cout << "2) Tridiagonal algorithm" << endl;
    cout << "3) Seidel and simple iterations methods" << endl;
    cout << "4) Method Jacobi" << endl;
    cout << "5) QR decomposition" << endl << endl;
    cout << "enter \"exit\" if you want to stop the programm" << endl << endl;

    int command = 8;
    while (command){
        cin >> command;
        if (command == 0)
        {
            cout << "Thank you for using this program. Come back again!";
            return 0;
        } else if (command == 1)
        {
            ifstream file_input("res/input_1.1.txt");
            int n;
            file_input >> n; 
            matrix coeffs(n, n), roots(n, 1);
            file_input >> coeffs >> roots;
            pair <matrix, matrix> lu_pair = lu(coeffs, roots, true);
            cout << "L =\n" << lu_pair.first << "\nU =\n" << lu_pair.second << endl;
            cout << "\nMatrix's determinant: " << det(coeffs, roots);
            cout << "\nEquations solving:\n" << solve(lu_pair, roots);
            cout << "\nInversed matrix:\n" << inverse(coeffs, roots) * coeffs;
            cout << endl << endl;
        }
        else if (command == 2)
        {
            cout << "Now this block in development" << endl;
        }  
        else if (command == 3)
        {
            cout << "Now this block in development" << endl;
        }
        else if (command == 4)
        {
            cout << "Now this block in development" << endl;
        }
        else if (command == 5)
        {
            cout << "Now this block in development" << endl;
        } else
            cout << "Sorry, but this command is invalid" << endl << endl;
    }
}
