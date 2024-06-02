#include <bits/stdc++.h>

using namespace std;
using double_matrix = vector<vector<double> >;


pair<int, int> get_indexes(const double_matrix& matrix1){
    double k = matrix1[0][1];
    pair<int, int> indexes = make_pair(0, 1);
    int n = matrix1.size(), m = matrix1[0].size();
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < m; j++)
            if (abs(matrix1[i][j]) > k) {
                k = abs(matrix1[i][j]);
                indexes = make_pair(i, j);
            }
    return indexes;
}


double ugol(double aji, double aii, double ajj){
    if (aii == ajj)
        return 3.14/4;
    return atan(2*aji/(aii-ajj)) / 2;
}


double_matrix mult(const double_matrix& matrix1, const double_matrix& matrix2) {
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


double_matrix ones_matr(int n){
    double_matrix E(n, vector<double>(n, 0));
    for(int i=0; i<n; i++)
        E[i][i] = 1;
    return E;
}


double_matrix u_matr(const double_matrix& matrix1){
    double_matrix res = ones_matr(matrix1.size());
    int i_max, j_max;
    tie(i_max, j_max) = get_indexes(matrix1);
    double phi = ugol(matrix1[i_max][j_max], matrix1[i_max][i_max], matrix1[j_max][j_max]);
    res[i_max][i_max] = cos(phi);
    res[j_max][j_max] = cos(phi);
    res[i_max][j_max] = -sin(phi);
    res[j_max][i_max] = sin(phi);
    return res;
}


double deviation(const double_matrix& matrix1){
    double eps = 0;
    int n=matrix1.size(), m=matrix1[0].size();
    for (int i=0; i<n; i++)
        for (int j=i+1; j<m; j++)
            eps += matrix1[i][j]*matrix1[i][j];
    return sqrt(eps);
}


double_matrix trans(const double_matrix& matrix1){
    int n = matrix1.size(), m = matrix1[0].size();
    double_matrix res(m, vector<double>(n));
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            res[j][i] = matrix1[i][j];
    return res;
}


pair<vector<double>, double_matrix> jacobi(double_matrix& coeff_matrix, double EPS){
    int n = coeff_matrix.size();
    double_matrix evect = ones_matr(n);
    while (EPS <= deviation(coeff_matrix)){
        double_matrix u = u_matr(coeff_matrix);
        evect = mult(evect, u);
        double_matrix trans_u = trans(u);
        coeff_matrix = mult(mult(trans_u, coeff_matrix), u);
    }

    vector<double> eval(n);
    for (int i=0; i<n; i++)
        eval[i] = coeff_matrix[i][i];
    return make_pair(eval, evect);
}


int main() {
    double_matrix coeff_matrix{{9, -2, 3}, {-2, 6, 8}, {3, 8, -6}};
    vector<double> eval;
    double_matrix evect;

    tie(eval, evect) = jacobi(coeff_matrix, 0.01);
    int n = eval.size();

    cout << "Собственные значения: " << endl;
    for (int i=0; i<n; i++)
        cout << "Значение " << i+1 << " = " << eval[i] << endl;
    cout << endl;

    cout << "Собственные векторы: " << endl;
    for (int i=0; i<n; i++){
        cout << "Вектор " << i+1 << " = ";
        for(double element: evect[i])
            cout << element << " ";
        cout << endl;
    }

    return 0;
}