#include "matrix.h"
#include <iostream>
#include <vector>



int main() {

    matrix A{
            {-3,1,-1},
            {6,9,-4},
            {5,-4,-8}


    };


    matrix res = solve(A,0.01);

    std::cout << "lambda :" << std::endl;
    print_matrix(res);


    return 0;
}
