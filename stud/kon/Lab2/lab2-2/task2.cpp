#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f1(double x1, double x2) {
    return x1-cos(x2)-3;
}

double f2(double x1, double x2){
    return x2-sin(x1)-3;
}

double df1dx1() {
    return 1;
}

double df1dx2(double x2) {
    return sin(x2);
}

double df2dx1(double x1) {
    return -cos(x1);
}

double df2dx2() {
    return 1;
}

double phi_1(double x2){
    return cos(x2)+3;
}

double phi_2(double x1){
    return sin(x1)+3;
}

pair<double, double> simple_iteration(double eps, double x1, double x2, int &iters) {
    while (iters < 1000) {
        double x1_new = phi_1(x2);
        double x2_new = phi_2(x1);
        if (fabs(x1_new - x1) < eps && fabs(x2_new - x2) < eps) {
            return {x1_new, x2_new};
        }
        x1 = x1_new;
        x2 = x2_new;
        iters++;
    }
    return {x1, x2};
}

pair<double, double> newton(double eps, double x1, double x2, int &iters) {
    while (iters < 1000) {
        double det = df1dx1() * df2dx2() - df1dx2(x2) * df2dx1(x1);
        double dx1 = (f1(x1, x2) * df2dx2() - f2(x1, x2) * df1dx2(x2)) / det;
        double dx2 = (f2(x1, x2) * df1dx1() - f1(x1, x2) * df2dx1(x1)) / det;
        x1 -= dx1;
        x2 -= dx2;
        if (fabs(dx1) < eps && fabs(dx2) < eps) {
            return {x1, x2};
        }
        iters++;
    }
    return {x1, x2};
}

int main(){
    double eps = 0.001;
    double x1 = 0.5, x2 = 0.5;
    int iterations = 0;
    ofstream fout("output.txt");
    auto solution_1 = simple_iteration(eps, x1, x2, iterations);
    fout << "Метод простых итераций:" << endl;
    fout << "x1 = " << solution_1.first << endl;
    fout << "x2 = " << solution_1.second << endl;
    fout << "Число итераций: " << iterations << endl;
    int iterations_2 = 0;
    auto solution_2 = newton(eps, x1, x2, iterations_2);
    fout << "Метод Ньютона:" << endl;
    fout << "x1 = " << solution_2.first << endl;
    fout << "x2 = " << solution_2.second << endl;
    fout << "Число итераций: " << iterations_2 << endl;
    return 0;
}