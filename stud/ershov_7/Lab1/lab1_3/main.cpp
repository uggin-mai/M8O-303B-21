#include "iteration.hpp"

using namespace std;

using matrix = matrix_t<double>;
using vec = vector<double>;

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix a(n);
    cin >> a;
    iter_solver my_solver(a, eps);
    vec b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vec x = my_solver.solve_simple(b);
    cout << "Simple iterations method" << endl;
    cout << "Iterations count: " << my_solver.iter_count << endl;
    cout << "System solution:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    vec z = my_solver.solve_zeidel(b);
    cout << "\nSeidel method" << endl;
    cout << "Iterations count: " << my_solver.iter_count << endl;
    cout << "System solution:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << z[i] << endl;
    }
}
