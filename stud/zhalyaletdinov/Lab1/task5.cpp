#include <bits/stdc++.h>

using namespace std;
using double_matrix = vector<vector<double>>;
using dvect = vector<double>;

double_matrix trans(const double_matrix& matrix1){
    int n = matrix1.size(), m = matrix1[0].size();
    double_matrix res(m, dvect(n));
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            res[j][i] = matrix1[i][j];
    return res;
}

double_matrix matr_plus(const double_matrix& matrix1, const double_matrix& matrix2) {
    double_matrix res;
    int n = matrix1.size();
    for (int i = 0; i < n; i++) {
        dvect row;
        for (int j = 0; j < n; j++) {
            row.push_back(matrix1[i][j] + matrix2[i][j]);
        }
        res.push_back(row);
    }

    return res;
}

double_matrix multiplone_matr(double_matrix& matrix1, double_matrix& matrix2) {
    int n1 = matrix1.size(), m1 = matrix1[0].size(), m2 = matrix2[0].size();
    double_matrix res(n1);
    for (int i=0; i<n1; i++)
        for (int j=0; j<m2; j++)
            res[i].push_back(0);

    for (int i=0; i<n1; i++) {
        for (int j=0; j<m2; j++) {
            double cntr = 0;
            for (int k=0; k<m1; k++)
               cntr += matrix1[i][k] * matrix2[k][j];
            res[i][j] = cntr;
        }
    }
    return res;
}

int sign(double a){
    return (a >= 0) ? 1 : -1;
}

double deviation(const double_matrix& matrix1){
    double eps = 0;
    int n = matrix1.size();
    for(int i=0; i<n; i++)
        for(int j=0; j<i-1; j++)
            eps += matrix1[i][j]*matrix1[i][j];
    return sqrt(eps);
}

double_matrix get_one_matr(int n){
    double_matrix E(n, dvect(n, 0));
    for(int i=0; i<n; i++)
        E[i][i] = 1;
    return E;
}

double_matrix get_H_matrix(const double_matrix& coeffts, int ind){
    int n = coeffts.size();
    double_matrix v(n, dvect(1));

    for(int i=0; i<n; i++){
        if (i < ind)
            v[i][0] = 0;
        else if (i == ind){
            double sum = 0;
            for (int j=ind; j < n; j++)
                sum += coeffts[j][i]*coeffts[j][i];
            v[i][0] = coeffts[i][i] + sign(coeffts[i][i]) * sqrt(sum);
        }
        else
            v[i][0] = coeffts[i][ind];   
    }

    double_matrix trans_v = trans(v);
    double k = -multiplone_matr(trans_v, v)[0][0]/2;
    v = multiplone_matr(v, trans_v);
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            v[i][j] /= k;

    double_matrix E = get_one_matr(n);
    return matr_plus(E, v);
}

pair<double_matrix, double_matrix> decompose(const double_matrix& coeff){
    double_matrix coeffts = coeff;
    double_matrix Q = get_H_matrix(coeffts, 0);
    coeffts = multiplone_matr(Q, coeffts);
    int n = coeffts.size();
    for (int i=1; i<n-1; i++){
        double_matrix H = get_H_matrix(coeffts, i);
        Q = multiplone_matr(Q, H);
        coeffts = multiplone_matr(H, coeffts);
    }
    return make_pair(Q, coeffts);
}

dvect evals_calc(double_matrix& coeffts, double EPS){
    while (deviation(coeffts) > EPS){
        pair<double_matrix, double_matrix> QR = decompose(coeffts);
        coeffts = multiplone_matr(QR.second, QR.first);
    }
    int n = coeffts.size();
    dvect r(n);
    for (int i=0; i<n; i++)
        r[i] = coeffts[i][i];
    return r;
}

int main() {
    double_matrix coefft_matrix{{-9, 2, 2}, {-2, 0, 7}, {8, 2, 0}};
    dvect evals = evals_calc(coefft_matrix, 0.01);
    for (int i=0; i<evals.size(); i++)
        cout << "Значение" << i+1 << ": " << evals[i] << endl;
}