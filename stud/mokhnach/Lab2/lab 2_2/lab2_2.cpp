#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
class matrix
{
    private:
        vector <vector <double>> obj;
    public:
        int  cols = 0, rows = 0;

        matrix() {}
        matrix(int _rows, int _cols)
        {
            rows = _rows;
            cols = _cols;
            obj = vector <vector <double>>(rows, vector <double>(cols));
        }

        operator double()
        {
            return obj[0][0];
        }

        vector <double>& operator[](int i)
        {
            return obj[i];
        }

        vector<double> getRow(int i){
            return obj[i];
        }
        void swap (int I, int J){
            for (int i = 0; i < cols; ++ i){
                swap (obj[I][i], obj[J][i]);
            }
        }

        matrix delete_column_row(int row,int coll) {
            int d_r = 0;
            matrix result = matrix (rows - 1, cols - 1);
            for (int i = 0; i < rows; ++i){
                if (i > row) d_r = 1;
                int d_c = 0;
                for (int j = 0; j < cols; ++j){
                    if (j > coll) d_c = 1;
                    if (i != row && j != coll){
                        result[i - d_r][j - d_c] = obj[i][j];
                    }
                }
            }
            return result;
        }

        double get_determinant() {
            double result = 0;
            if (rows == 1){ 
                return obj[0][0];  
            }
            for (int k = 0; k < rows; ++k) {
                matrix M = matrix(rows - 1, rows - 1);
                for (int i = 1; i < rows; ++i) {
                    int t = 0;
                    for (int j = 0; j < rows; ++j) {
                    if (j == k)
                        continue;
                    M[i-1][t] = obj[i][j];
                    t += 1;
                    }
                }
                result += pow(-1, k + 2) * obj[0][k] * M.get_determinant();
            }
            return result;
        }
        double get_norm() {
            double result = 0;
            for (int i = 0; i < rows; ++i){
                for (int j = 0; j < cols; ++j){
                    result += pow(obj[i][j], 2);
                }
            }
            return pow(result, 0.5);
        }
};



//****Определение операторов****//
istream& operator>>(istream& stream, matrix& m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++)
            stream >> m[i][j];
    }
    return stream;
}

ostream& operator<<(ostream& stream, matrix m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++)
            stream << m[i][j] << ' ';
        stream << '\n';
    }
    return stream;
}

matrix operator*(double a, matrix b)
{
    for (int i = 0; i < b.rows; i++)
    {
        for (int j = 0; j < b.cols; j++)
            b[i][j] *= a;
    }
    return b;
}

matrix operator+(matrix a, matrix b)
{
    if (a.rows != b.rows || b.cols != a.cols)
        return matrix(0, 0);
    matrix res(a.rows, a.cols);
    for (int i = 0; i < b.rows; i++)
    {
        for (int j = 0; j < res.cols; j++)
            res[i][j] = a[i][j] + b[i][j];
    }
    return res;
}

matrix operator-(matrix a, matrix b)
{
    if (a.rows != b.rows || b.cols != a.cols)
        return matrix(0, 0);
    matrix res(a.rows, a.cols);
    for (int i = 0; i < b.rows; i++)
    {
        for (int j = 0; j < res.cols; j++)
            res[i][j] = a[i][j] - b[i][j];
    }
    return res;
}

matrix operator*(matrix a, matrix b)
{
    if (a.cols != b.rows)
        return matrix(0, 0);
    matrix res(a.rows, b.cols);
    for (int i = 0; i < res.rows; i++) 
    {
        for (int j = 0; j < res.cols; j++)
        {
            res[i][j] = 0;
            for (int k = 0; k < a.cols; k++)
                res[i][j] += a[i][k] * b[k][j];
        }
    }
    return res;
}

matrix get_J(matrix X){
    matrix J = matrix(2, 2);
    J[0][0] = X[0][0]/8;
    J[0][1] = X[1][0]/2;
    J[1][0] = -exp(X[0][0]) - 1;
    J[1][1] = 4;
    return J;
}

matrix get_A1(matrix X){
    matrix A = matrix(2, 2);
    A[0][0] = pow(X[0][0], 2)/16 + pow(X[1][0], 2)/4 - 1;
    A[0][1] = X[1][0]/2;
    A[1][0] = 4*X[1][0] - exp(X[0][0]) -  X[0][0];
    A[1][1] = 4;
    return A;
}

matrix get_A2(matrix X){
    matrix A = matrix(2, 2);
    A[0][0] = X[0][0]/8;
    A[0][1] = pow(X[0][0], 2)/16 + pow(X[1][0], 2)/4 - 1;
    A[1][0] = -exp(X[0][0]) - 1;
    A[1][1] = 4*X[1][0] - exp(X[0][0]) -  X[0][0];
    return A;
}

matrix Newton_method(matrix X0, double eps, int &iterations) {
    matrix X1 = matrix(2, 1);
    X1[0][0] = 1; X1[1][0] = 1;
    matrix dets_div_result = matrix(2, 1);
    while ((X1 - X0).get_norm() > eps){
        X0 = X1;
        dets_div_result[0][0] = get_A1(X0).get_determinant()/get_J(X0).get_determinant();
        dets_div_result[1][0] = get_A2(X0).get_determinant()/ get_J(X0).get_determinant();
        X1 = X0 - dets_div_result;
        iterations += 1;
    }
    return X1;
}

matrix Phi(matrix X){
    matrix phi = matrix(2,1);
    phi[0][0] = 4*X[1][0] - exp(X[0][0]);
    phi[1][0] =  (1 - pow(X[0][0], 2)/16)*4/X[1][0];
    // phi[1][0] = (16 - X[1][0]*X[1][0]*4)/X[0][0];
    // phi[0][0] = (exp(X[0][0]) + X[0][0])/4;
    return phi;
}

matrix Iterations_method(matrix X0, double eps, int &iterations) {
    matrix X1 = matrix(2, 1);
    X1[0][0] = 1.709; X1[1][0] = 1.808;
    double q = 0.4;
    matrix phi = matrix(2, 1);
    matrix diff = X1 - X0;
    do{
        X0 = X1;
        X1 = Phi(X0);
        iterations += 1;
        diff = X1 - X0;
        if (iterations > 1000) break;
    } while (q / (1 - q) * diff.get_norm() > eps);
    return X1;
}

int main() {
    cout.precision(9);
    ofstream fout("answer.txt");
    ifstream fin("input.txt");
    double eps;
    fin >> eps;
    matrix X0_iterations(2,1);
    cin >> X0_iterations[0][0] >> X0_iterations[1][0];
    int iterations_iterations = 0;
    matrix X_iterations = Iterations_method(X0_iterations, eps, iterations_iterations);
    cout << "===Iterations method===\nIterations number: " << iterations_iterations << "\nRoots:\n" << X_iterations << '\n';

    int iterations_Newton = 0;
    matrix X0_Newton(2,1);
    X0_Newton[0][0] = 1.5; X0_Newton[1][0] = 1.5;
    matrix X_Newton = Newton_method(X0_Newton, eps, iterations_Newton);
    fout << "===Newton method===\nIterations number: " << iterations_Newton << "\nRoots:\n" << X_Newton << '\n';
}