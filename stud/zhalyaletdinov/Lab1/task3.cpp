#include <bits/stdc++.h>

using namespace std;
using double_matrix = vector<vector<double>>;


double_matrix matr_plus(const double_matrix& matrix1, const double_matrix& matrix2) {
    int n = matrix1.size();
    int m = matrix1[0].size();

    double_matrix res(n, vector<double>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    return res;
}


double_matrix mult(const double_matrix& matrix1, const double_matrix& matrix2) {
    int n1 = matrix1.size(), m1 = matrix1[0].size(), n2 = matrix2.size(), m2 = matrix2[0].size();
    double_matrix res(n1, vector<double>(m2, 0));

    for (int i = 0; i < n1; ++i)
        for (int j = 0; j < m2; ++j) {
            double cntr = 0;
            for (int k = 0; k < m1; ++k)
                cntr += matrix1[i][k] * matrix2[k][j];
            res[i][j] = cntr;
        }

    return res;
}


void matrix_preprocess(double_matrix& coeffs, double_matrix& r) {
    double n = coeffs.size(), m = coeffs[0].size();
    for(int i=0; i < n; i++) {
        double a = coeffs[i][i];
        if (a == 0)
            continue;
        r[i][0] /= a;
        for (int j=0; j < m; j++)
            if (i == j)
                coeffs[i][j] = 0;
            else
                coeffs[i][j] /= -a;
    }
}


double eps_curr(const double_matrix& vect1, const double_matrix& vect2) {
    double eps = 0;
    for(int i=0; i<vect1.size(); i++)
        eps += pow(vect1[i][0] - vect2[i][0], 2);
    return sqrt(eps);
}


void cout_matrix(const double_matrix& matrix1) {
    for(const auto& vect: matrix1) {
        cout << '\t';
        for (auto x: vect)
            cout << x << " ";
        cout << endl;
    }
}


pair<int, double_matrix> iterations(double_matrix& a, double_matrix& b, double EPS0) {
    double_matrix p_x = b, c_x;
    int k = 0;
    while (eps_curr(c_x = matr_plus(b, mult(a, p_x)), p_x) > EPS0) {
        k += 1;
        p_x = c_x;
    }
    return {k, p_x};
}


pair<int, double_matrix> seidel(const double_matrix& a, const double_matrix& b, double EPS0) {
    double_matrix p_x = b, c_x = p_x;
    bool flag = true;
    int k = 0;
    double eps = 0;
    while (flag or eps > EPS0) {
        k += 1;
        flag = false;
        for(int i = 0; i < b.size(); i++) {
            c_x[i][0] = 0;
            for(int j=0; j < b.size(); j++) {
                if (i == j)
                    c_x[i][0] += b[i][0];
                else if (i < j)
                    c_x[i][0] += a[i][j]*p_x[j][0];
                else
                    c_x[i][0] += a[i][j]*c_x[j][0];
            }
        }
        eps = eps_curr(c_x, p_x);
        p_x = c_x;
    }
    return {k, p_x};
}


int main() {
    double_matrix coeff_matrix = {
        {-7, -1, 2, 2},
        {3, -20, 0, -8},
        {-9, 1, 18, -6},
        {-1, 0, -1, -6}
    };

    double_matrix roots = {
        {-24},
        {-47},
        {28},
        {-50}
    };

    matrix_preprocess(coeff_matrix, roots);

    int last_iter;
    double_matrix res_matr;

    tie(last_iter, res_matr) = iterations(coeff_matrix, roots, 0.001);
    cout << "Метод простых итераций" << endl << "\tКоличество итераций: " << last_iter << endl << "\tРешения СЛАУ: " << endl;
    cout_matrix(res_matr);
    
    cout << endl << endl;

    tie(last_iter, res_matr) = seidel(coeff_matrix, roots, 0.001);
    cout << "Метод Зейделя" << endl << "\tКоличество итераций: " << last_iter << endl << "\tРешения СЛАУ: " << endl;
    cout_matrix(res_matr);
}
