#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;



class Matrix {
	public:
	int rows, cols;
   double **a;

   Matrix (){
      rows = 0;
      cols = 0;
      a = new double*[rows];
      for(int i = 0; i < rows; i++)
         a[i] = new double[cols];
   }
    
    Matrix(int n, int m, bool identity = 0)
	{
      rows = n;
      cols = m;
		a = new double *[rows];
      for (int i = 0; i < rows; i++){
         a[i] = new double[cols];
         for (int j = 0; j < cols; j++){
            a[i][j] = identity * (i == j);
         }
      }
	}

    double* operator[](int i)
	{
      return a[i];
   	}

    double norm(){
        double norm = 0;
        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                norm += pow(a[i][j], 2);
            }
        }
        return pow(norm, 0.5);
    }

    void print()
    {
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                cout <<"\t"<< a[i][j] << "\t";
            }
            cout << endl;
        }
    }
};

Matrix operator* (Matrix & A, double k){
   Matrix C = Matrix (A.rows, A.cols);
    for (int i = 0; i < A.rows; ++ i){
        for (int j = 0; j < A.cols; ++ i){
            C[i][j] = A[i][j] * k;
        }
    }
    return C;
}

Matrix operator* (Matrix & A, Matrix & B){
   Matrix C = Matrix (A.rows, B.cols);
   for (int i = 0; i < A.rows; ++ i){
      for (int j = 0; j < B.cols; ++ j){
         for (int k = 0; k < A.cols; ++ k){
            C[i][j] += A[i][k] * B[k][j];
         }
      }
   }
   return C;
}

Matrix operator+ (Matrix A, Matrix B){
   Matrix C = Matrix(A.rows, A.cols);
   for (int i = 0; i < A.rows; ++i){
      for (int j = 0; j < A.cols; ++j){
         C[i][j] = A[i][j] + B[i][j];
      }
   }
   return C;
}

Matrix operator- (Matrix A, Matrix B){
   Matrix C = Matrix(A.rows, A.cols);
   for (int i = 0; i < A.rows; ++i){
      for (int j = 0; j < A.cols; ++j){
         C[i][j] = A[i][j] - B[i][j];
      }
   }
   return C;
}

double det(Matrix A){
    double d = 0;
    int n = A.rows;
    if (n == 1){ 
        return A[0][0];  
    }
    else if (n == 2){
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }
    else {
        for (int k = 0; k < n; ++k) {
            Matrix M = Matrix(n - 1, n - 1);
            for (int i = 1; i < n; ++i) {
                int t = 0;
                for (int j = 0; j < n; ++j) {
                if (j == k)
                    continue;
                M[i-1][t] = A[i][j];
                t += 1;
                }
            }
            d += pow(-1, k + 2) * A[0][k] * det(M);
        }
    }
    return d;
}

double Norm(Matrix A){
    double norm = 0;
    for (int i = 0; i < A.rows; ++i){
        for (int j = 0; j < A.cols; ++j){
            norm += pow(A[i][j], 2);
        }
    }
    return pow(norm, 0.5);
}

Matrix Minor(Matrix A, int a, int b){
    int n = A.rows;
    int m = A.cols;
    Matrix C = Matrix(n - 1, m - 1);
    int k = 0;
    for (int i = 0; i < n; ++i){
        if (i > a){
            k = 1;
        }
        for (int j = 0; j < m; ++j){
            if (i != a){
                if (j < b){
                    C[i - k][j] = A[i][j];
                }
                else if (j > b){
                    C[i - k][j - 1] = A[i][j];
                }
            }
        }
    }
    return C;
}

Matrix Transpose(Matrix A){
   int n = A.rows;
   int m = A.cols;
   Matrix C = Matrix (m, n);
   for (int i = 0; i < n; ++ i){
      for (int j = 0; j < m; ++ j){
         C[j][i] = A[i][j];
      }
   }
   return C;
} 

Matrix Inverse(Matrix A){
   double d = det(A);
   int n = A.rows;
   Matrix inv_A = Matrix(n, n);
   for(int i = 0; i < n; ++i){
      for(int j = 0; j < n; ++j){
         Matrix M = Matrix(n - 1, n - 1);
         M = Minor(A, i, j);
         inv_A[i][j] = pow(-1.0, i + j + 2) * det(M) / d;
      }
   }
   inv_A = Transpose(inv_A);
   return inv_A;
}

Matrix Jacobian(Matrix X, int k = -1)
    {
        Matrix J = Matrix(2, 2);
        J[0][0] = 2 * X[0][0] * X[1][0];
        J[0][1] = X[0][0] * X[0][0] + 4;
        J[1][0] = 2 * X[0][0] - 2;
        J[1][1] = 2 * X[1][0] - 2;
        if (k != -1){
            J[k][0] = (pow(X[0][0], 2) + 4) * X[1][0] - 8;
            J[k][1] = pow((X[0][0] - 1), 2) + pow((X[1][0] - 1), 2) - 4;
        }
        return J;
    }

Matrix F(Matrix X){
    Matrix f = Matrix(2,1);
    f[0][0] = (pow(X[0][0], 2)+4)*X[1][0]-8;
    f[1][0] = pow(X[0][0]-1, 2)+pow(X[1][0]-1, 2)-4;
    return f;
}

Matrix EqF(Matrix X){
    Matrix eqf = Matrix(2,1);
    //double x1 = X[0][0];
    //double x2 = X[1][0];
    eqf[0][0] = 8 / (pow(X[0][0], 2) + 4);
    eqf[1][0] = pow(4 - pow(X[1][0]-1, 2), 0.5) + 1;
    return eqf;
}

int main()
{
    // Ньютон
    double eps = 0.001;
    Matrix X[2] = {Matrix(2, 1), Matrix(2, 1)};
    X[0][0][0] = 2.8; X[0][1][0] = 0.6;
    X[1][0][0] = 3.0; X[1][1][0] = 0.8;
    int k = 0;
    Matrix J = Matrix(2, 1);
    Matrix f = Matrix(2, 1);
    while (Norm(X[1] - X[0]) > eps){
        X[0] = X[1];
        J = Inverse(Jacobian(X[1]));
        f = F(X[1]);
        X[1] = X[1] - (J * f);
        k += 1;
    }
    X[1].print();
    cout << "\n" << k << "\n";

    // Простые итерации
    k = 0;
    X[0][0][0] = 2.8; X[0][1][0] = 0.6;
    X[1][0][0] = 3.0; X[1][1][0] = 0.8;
    double q = 0.4;
    Matrix eqf = Matrix(2, 1);
    while (q / (1 - q) * Norm(X[1] - X[0]) > eps and k < 10){
        X[0] = X[1];
        X[1] = EqF(X[1]);
        //X[1].print();
        //cout << "\n";
        k += 1;
    }
    X[1].print();
    cout << k << "\n";

}
