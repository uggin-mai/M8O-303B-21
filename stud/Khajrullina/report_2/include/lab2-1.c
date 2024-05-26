#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const double a = 0;

double Function(double x) {
    return log10(x + 1) - x + 0.5;
}

double Derivative(double x) {
    return (1 / (log(10) * (x + 1))) - 1;
}

double simple_iteration_function(double x) {
    return log10(x + 1) + 0.5;
}

double absolute(double a) {
    return a > 0 ? a : -a;
}

double simple_iteration_method(double (*Function)(double), double eps) {
    const double constanta = 0.5;
    double result, prev = a;
    int iter;
    for (iter = 1; eps < constanta / (1 - constanta) * absolute((result = Function(prev)) - prev); iter++) {
        prev = result;
    }
    printf("Iterations: %d\n", iter);
    return result;
}

double newton_method(double (*Function)(double), double (*Derivative)(double), double eps) {
    double result, prev = a;
    int iter;
    for (iter = 1; eps < absolute((result = prev - Function(prev) / Derivative(prev)) - prev); iter++) {
        prev = result;
    }
    printf("Iterations: %d\n", iter);
    return result;
}

int main(void) {
    float eps;

    scanf("%f", &eps);

    if (eps <= 0) {
        fprintf(stderr, "Error\n");
        return 0;
    }
    printf("Newton method:\n");
    printf("x* = %.4f\n", newton_method(Function, Derivative, eps));
    printf("Simple iteration method:\n");
    printf("x* = %.4f\n", simple_iteration_method(simple_iteration_function, eps));

    return 0;
}