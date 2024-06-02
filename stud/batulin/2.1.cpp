#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

double initfunc(double x){
    return pow(2, x) - pow(x, 2) - 0.5;
}

double dif(double x){
    return pow(2, x) * log(2) - (2 * x);
}

double ddif(double x){
    return pow(log(2),2) * pow(2, x) - 2;
}

double eqf(double x){
    return pow(pow(2, x) - 0.5, 0.5);
}

double deqf(double x){
    return log(2) * pow(2, x - 1/2) / sqrt(pow(2,x + 1) - 1);
}

int main(){

    // Ньютон
    double x0 = 1.6;
    double eps = 0.001;
    if (initfunc(x0) * ddif(x0) > 0){
        cout << "Начальная точка подходит" << "\n";

        int k = 0;
        double x[2] = {x0 - eps - 1, x0};
        while (abs(x[1] - x[0]) >= eps){
            x[0] = x[1];
            x[1] -= initfunc(x[1]) / dif(x[1]);
            k += 1;
        }
        cout << x[1] << " " << k << "\n";
    }
    else{
        cout << "Начальная точка не подходит" << "\n";
    }

    // Простые итерации
    int k = 0;
    double q = abs(deqf(1.6));
    double x[2] = {1.5, (1.5 + 1.6) / 2};
    while (q / (1 - q) * abs(x[1] - x[0]) > eps){
        x[0] = x[1];
        x[1] = eqf(x[1]);
        k += 1;
    }
    cout << x[1] << " " << k << "\n";

}
