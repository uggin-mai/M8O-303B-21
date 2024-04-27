#include <bits/stdc++.h>

using namespace std;

pair<int, double> iteration_solution(double x_start, double q, double e){
    double r = x_start, l = x_start*5;
    int n=0;

    while (q*abs(r - l)/(1-q) >= e){
        n++;
        l = r;
        r = pow(log(r+1) + 1, 1.0/3);
    }
    return {n, r};
}

pair<int, double> newton_solution(double x_start, double e){
    double r = x_start - (log(x_start+1) - pow(x_start, 3) + 1)/(1/(x_start+1) - 3*pow(x_start, 2)), l = x_start;
    int n=0;

    while (abs(r - l) >= e){
        n++;
        l = r;
        r = r - (log(r+1) - pow(r, 3) + 1)/(1/(r+1) - 3*pow(r, 2));
    }
    return {n, r};
}


int main(){
    pair<int, double> newton_algo = newton_solution(0.5, 0.000001), iteration_algo = iteration_solution(0.5, 0.5, 0.000001);
    cout << endl << "log(x+1) - x^3 + 1 = 0" << endl << endl << "k = " << newton_algo.first << endl << "x = " << newton_algo.second << endl << endl << "k = " << iteration_algo.first << endl << "x = " << iteration_algo.second << endl << endl;
    return 1;
}
