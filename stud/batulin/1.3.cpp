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

	void swapRows(int row1, int row2)
	{
		swap(a[row1], a[row2]);
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

   Matrix minor(int y, int z)
   {
      Matrix C = Matrix (rows - 1, cols - 1);
      int k = 0;
      for (int i = 0; i < rows; ++i){
         if (i > y){
            k = 1;
         }
         for (int j = 0; j < cols; ++j){
            if (i != y){
               if (j < z){
                  C[i - k][j] = a[i][j];
               }
               else if (j > z){
                  C[i - k][j - 1] = a[i][j];
               }
            }
         }
      }
      return C;
   }

   Matrix inverse()
   {
      double d = this->determinant();
      Matrix inv_A = Matrix(rows, cols);
      for(int i = 0; i < rows; ++i){
         for(int j = 0; j < rows; ++j){
            Matrix M = Matrix(rows - 1, rows - 1);
            M = this->minor(i, j);
            inv_A[i][j] = pow(-1.0, i + j + 2) * M.determinant() / d;
         }
      }
      inv_A = inv_A.transpose();
      return inv_A;
   }

   double determinant()
   {
      double d = 0;
      if (rows == 1){ 
         return a[0][0];  
      }
      else if (rows == 2){
         return a[0][0] * a[1][1] - a[0][1] * a[1][0];
      }
      for (int k = 0; k < rows; ++k) {
            Matrix M = Matrix(rows - 1, cols - 1);
            for (int i = 1; i < rows; ++i) {
               int t = 0;
               for (int j = 0; j < rows; ++j) {
                  if (j == k)
                     continue;
                  M[i-1][t] = a[i][j];
                  t += 1;
               }
            }
            d += pow(-1, k + 2) * a[0][k] * M.determinant();
         }
      return d;
   }
   double Norm(){
      double norm = 0;
      for (int i = 0; i < rows; ++i){
         for (int j = 0; j < cols; ++j){
            norm += pow(a[i][j], 2);
         }
      }
      return pow(norm, 0.5);
   }
};

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

/*
19 -4 -9 -1 -2 20 -2 -7 6 5 -25 9 0 -3 -9 12
100 -5 34 69
*/

int main()
{

   int n = 4;
   Matrix A = Matrix (n, n);
   Matrix Alpha = Matrix (n, n);
   Matrix Down = Matrix (n, n);
   Matrix E = Matrix (n, n, 1);
   Matrix B = Matrix (n, 1);
   Matrix Beta = Matrix(n, 1);
   A.scan();
   B.scan();

   // приведение СЛАУ у эквивалентному виду
   for (int i = 0; i < n; ++i){
      for (int j = 0; j < n; ++j){
         double s = 0;
         if (i != j){
            Alpha[i][j] = -A[i][j] / A[i][i];
            if (i > j){
               Down[i][j] = Alpha[i][j];
            }
            s += abs(A[i][j]);
         }
         else{
            Alpha[i][j] = 0;
         }
      }
      Beta[i][0] = B[i][0] / A[i][i];
   }

   // решение СЛАУ
   Matrix X[2] = {Beta, Beta + Alpha * Beta};
   int k = 1;
   double eps = 0.01;
   while ((X[1] - X[0]).Norm() > eps){
      X[0] = X[1];
      X[1] = Beta + Alpha * X[1];
      k += 1;
   }
   X[1].print();
   cout << "\n";
   cout << k << "\n";

   // решение СЛАУ методом Зейделя
   k = 1;
   Matrix C = E - Down;
   C = C.inverse();
   C = C * Beta;
   X[0] = Beta;
   X[1] = Alpha * Beta + C;
   while ((X[1] - X[0]).Norm() > eps){
      X[0] = X[1];
      X[1] = Alpha * X[1] + C;
      k += 1;
   }

   cout << "\n";
   X[1].print();
   cout << k << "\n";
}
