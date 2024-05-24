#include <bits/stdc++.h>

using namespace std;
using double_matrix = vector<vector<double>>;


double_matrix tridiag(double_matrix& coeffs, double_matrix& res_matr) {
    double a, b, c, d;
    a = 0;
    b = coeffs[0][0];
    c = coeffs[0][1];
    d = res_matr[0][0];
    vector<double> P(coeffs[0].size(), 0), Q(coeffs[0].size(), 0);

    P[0] = -c/b;
    Q[0] = d/b;
    for (int i=1; i < coeffs.size() - 1; i++){
        a = coeffs[i][i-1];
        b = coeffs[i][i];
        c = coeffs[i][i+1];
        d = res_matr[i][0];

        P[i] = -c/(b + a*P[i-1]);
        Q[i] = (d - a*Q[i-1])/(b + a*P[i-1]);
    }

    a = coeffs[coeffs.size()-1][coeffs[0].size()-2];
    b = coeffs[coeffs.size()-1][coeffs[0].size()-1];
    c = 0;
    d = res_matr[res_matr.size()-1][0];

    Q[Q.size()-1] = (d - a * Q[Q.size()-2]) / (b + a * P[P.size()-2]);

    double_matrix decision(res_matr.size());
    for(int i=0; i<decision.size(); i++)
        decision[i].push_back(0);

    decision[decision.size()-1][0] = Q[Q.size()-1];
    for (int i = decision.size()-2; i > -1; i--)
        decision[i][0] = P[i]*decision[i+1][0] + Q[i];

    return decision;
}


void cout_matrix(double_matrix& matrix1) {
    for(const auto& vect: matrix1) {
        for (auto x: vect)
            cout << x << " ";
        cout << endl;
    }
}


int main() {
    double_matrix coeff_matrix{{-11, -8, 0, 0, 0}, {9, -17, 1, 0, 0}, {0, -4, 20, 9, 0}, {0, 0, -4, -14, 3}, {0, 0, 0, -6, 14}};
    double_matrix decisions = {{99}, {-75}, {66}, {54}, {8}};
    
    double_matrix res_matrix = tridiag(coeff_matrix, decisions);
    cout << "Решения СЛАУ: " << endl;
    for (int i=0; i<res_matrix.size(); i++)
        printf("x%d = %f \n", i+1, res_matrix[i][0]);

    return 0;
}