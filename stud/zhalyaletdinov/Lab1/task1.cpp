#include <bits/stdc++.h>

using namespace std;
using double_matrix = vector<vector<double> >;


pair<double_matrix, double_matrix> decompose(double_matrix& coeffs, double_matrix& roots) {
    int shape=coeffs.size(), shape2=coeffs[0].size(), shape3 = roots[0].size();
    double_matrix L(shape), U = coeffs;
    for (int i=0; i<shape; i++)
        for (int j=0; j<shape2; j++)
            L[i].push_back(0);

    for (int k=0; k<shape; k++) {
        if (U[k][k] == 0) {
            for (int i=k+1; i<shape; i++) {
                if (U[i][k] != 0) {
                    swap(U[k], U[i]);
                    swap(L[k], L[i]);
                    swap(coeffs[k], coeffs[i]);
                    swap(roots[k], roots[i]);
                    break;
                }
            }
        }
        L[k][k] = 1;
        for (int i=k+1; i<shape; i++) {
            L[i][k] = U[i][k]/U[k][k];
            if (U[i][k] == 0)
                continue;
            for(int j=k; j<shape2; j++)
                U[i][j] -= L[i][k]*U[k][j];

        }
    }

    return {L, U};
}


double operdelitel(const double_matrix& U) {
    double det = 1;
    for (int i=0; i<U.size(); i++)
        det *= U[i][i];
    return det;
}


double_matrix calculate_decisions(double_matrix& coeffs, double_matrix& roots) {
    auto [L, U] = decompose(coeffs, roots);
    double_matrix res = roots;

    for (int k=0; k<res[0].size(); k++)
        for (int i=0; i<res.size(); i++)
            for (int j=0; j<i; j++)
                res[i][k] -= res[j][k]*L[i][j];
    for (int k=0; k<res[0].size(); k++) {
        for (int i=coeffs.size()-1; i>-1; i--) {
            for (int j=i+1; j<roots.size(); j++) {
                res[i][k] -= res[j][k]*U[i][j];
            }
            res[i][k] /= U[i][i];
        }
    }

    return res;
}


double_matrix get_inverse_matrix(double_matrix& matrix1) {
    double_matrix E(matrix1.size());
    for (int i=0; i<matrix1.size(); i++)
        for (int j=0; j<matrix1.size(); j++)
            E[i].push_back((i == j) ? 1 : 0);
    return calculate_decisions(matrix1, E);
}

void print_matrix(const double_matrix& matrix1) {
    for(const auto& vect: matrix1) {
        for (auto val: vect)
            cout << val << " ";
        cout << endl;
    }
}

int main() {
    double_matrix coefficient_matrix{{-4, -9, 4, 3}, {2, 7, 9, 8}, {4, -4, 0, -2}, {-8, 5, 2, 9}};
    double_matrix equation_roots = {{-51}, {76}, {26}, {-73}};

    auto [l, u] = decompose(coefficient_matrix, equation_roots);

    cout << endl << "L матрица: " << endl;
    print_matrix(l);

    cout << endl << endl;
    cout << "U матрица: " << endl;
    print_matrix(u);

    cout << endl << endl;
    cout << "Определитель матрицы равен: " << operdelitel(u) << endl;

    cout << endl << endl << "Обратная матрица: " << endl;
    double_matrix inversed = get_inverse_matrix(coefficient_matrix);
    print_matrix(inversed);

    cout << endl << endl << "Решения системы: " << endl;
    double_matrix decisions = calculate_decisions(coefficient_matrix, equation_roots);
    print_matrix(decisions);

    return 0;
}