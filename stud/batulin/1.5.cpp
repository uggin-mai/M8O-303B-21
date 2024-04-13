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

   void scan()
   {
      for(int i = 0; i < rows; i++)
      {
         for(int j = 0; j < cols; j++)
         {
            scanf ("%lf", &a[i][j]);
         }
      }
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

	double* operator[] (int i)
	{
      return a[i];
   }
   
   void swapRows(int row1, int row2)
	{
		swap(a[row1], a[row2]);
	}

   Matrix transpose()
   {
      Matrix C = Matrix (cols, rows);
      for (int i = 0; i < rows; ++ i){
         for (int j = 0; j < cols; ++ j){
            C[j][i] = a[i][j];
         }
      }
      return C;
   }

   double crit()
   {
      double s = 0;
      for (int i = 0; i < rows - 1; ++i){
         for (int j = i + 1; j < cols; ++j){
            s += pow(a[i][j], 2);
         }
      }
      return pow(s, 0.5);
   }

   pair<int, int> find_max(){
   int i_max = 0;
   int j_max = 1;
   for (int i = 0; i < rows - 1; ++i){
      for (int j = i + 1; j < rows; ++j){
         if (abs(a[i][j]) > abs(a[i_max][j_max])){
            i_max = i;
            j_max = j;
         }
      }
   }
   return pair<int, int>(i_max, j_max);
}
};

Matrix operator* (Matrix & A, double mult){
   Matrix C = Matrix (A.rows, A.cols);
   for (int i = 0; i < A.rows; ++i)
   {
      for (int j = 0; j < A.cols; ++j)
      {
         C[i][j] += A[i][j] * mult;
      }
   }
   return C;
}

Matrix operator* (Matrix & A, Matrix & B){
   Matrix C = Matrix (A.rows, B.cols);
   for (int i = 0; i < A.rows; ++i){
      for (int j = 0; j < B.cols; ++j){
         for (int k = 0; k < A.cols; ++k){
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

double scalar(Matrix &A, Matrix &B){
    double s;
    for (int i = 0; i < A.rows; ++i){
        for (int j = 0; j < B.cols; ++j){
            for (int k = 0; k < A.cols; ++k){
                s += A[i][k] * B[k][j];
            }
        }
    }
    return s;
}

int sign(double a){
    return (a >= 0) ? 1 : -1;
}

pair<Matrix, Matrix> QR_decomposition(Matrix A){
    int n = A.rows;
    Matrix E = Matrix(n, n, 1);
    Matrix Q = Matrix(n, n, 1);
    Matrix H = Matrix(n, n);
    for (int j = 0; j < n - 1; ++j){
        Matrix v = Matrix(n, 1);
        for (int i = j; i < n; ++i){
            v[i][0] = A[i][j];
            if (i == j){
                double s = 0;
                for (int k = j; k < n; ++k){
                    s += pow(A[k][j], 2);
                }
                v[i][0] += sign(A[i][j]) * pow(s, 0.5);
            }
        }
        Matrix v_T = v.transpose();
        Matrix v1 = v * v_T;
        H = v1 * 2;
        H = E - H * (1 / scalar(v_T, v));
        Q = Q * H;
        A = H * A;
    }
    return pair<Matrix, Matrix> (Q, A);
}

pair<double, double> complex(double b, double c){
    if (pow(b, 2) - 4 * c < 0){
        return pair<double, double> (-b / 2, pow(- pow(b, 2) + 4 * c, 0.5) / 2);
    }
    else{
        return pair<double, double> (-b / 2 + pow(pow(b, 2) - 4 * c, 0.5) / 2, 0);
    }
}

bool endcheck(Matrix A[2], double eps){
    int n = A[1].rows;
    int j = 0;
    while (j < n - 1){
        pair<double, double> z0 = complex( - A[0][j][j] - A[0][j + 1][j + 1], A[0][j][j] * A[0][j + 1][j + 1] - A[0][j][j + 1] * A[0][j + 1][j]);
        pair<double, double> z1 = complex( - A[1][j][j] - A[1][j + 1][j + 1], A[1][j][j] * A[1][j + 1][j + 1] - A[1][j][j + 1] * A[1][j + 1][j]);
        if (z1.second == 0){
            double s = 0;
            for (int i = j + 1; i < n; ++i){
                s += pow(A[1][i][j], 2);
            }
            s = pow(s, 0.5);
            if (s > eps){
                return true;
            }
        }
        else{
            if (pow(pow(z1.first - z0.first, 2) + pow(z1.second - z0.second, 2), 0.5) > eps){
                return true;
            }
            else{
                j += 1;
            }
        }
        j += 1;
    }
    return false;
}

/*
3 -7 -1 -9 -8 7 5 2 2
*/

int main()
{
    int n = 3;
    Matrix A = Matrix (n, n);
    A.scan();

    // начало QR-разложения
    double eps = 0.01;
    Matrix AA[2] = {A, Matrix(n, n)};
    auto [Q, R] = QR_decomposition(A);
    AA[1] = R * Q;
    int k = 1;

    while (endcheck(AA, eps)){
        auto [Q, R] = QR_decomposition(AA[1]);
        AA[0] = AA[1];
        AA[1] = R * Q;
        k += 1;
    }

    // получение собственных значений
    string eigenvalue[n];
    int i = 0;
    while (i < n - 1){
        pair<double, double> z = complex(AA[1][i][i] + AA[1][i + 1][i + 1], AA[1][i][i] * AA[1][i + 1][i + 1] - AA[1][i][i + 1] * AA[1][i + 1][i]);
        if (z.second == 0){
            eigenvalue[i] = to_string(AA[1][i][i]);
        }
        else{
            eigenvalue[i] = to_string(z.first) + " + " + to_string(z.second) + "i";
            i += 1;
            eigenvalue[i] = to_string(z.first) + " - " + to_string(z.second) + "i";
        i += 1;
        }
        i += 1;
    }

    AA[1].print();
    cout << k << "\n";
    for (int i = 0; i < n; ++i){
        cout << eigenvalue[i] << " ";
    }
}
