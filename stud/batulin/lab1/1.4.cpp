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
-7 4 5 4 -6 -9 5 -9 -8
*/

int main()
{
   int n = 3;
   Matrix A = Matrix(n, n);
   A.scan();
   Matrix X = Matrix(n, n, 1);
   double eps = 0.01;
   int k = 0;

   // метод вращений Якоби
   while (A.crit() > eps){
      Matrix U = Matrix(n, n, 1);
      Matrix U_T = Matrix(n, n, 1);
      auto [i_max, j_max] = A.find_max();
      double phi = atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max])) / 2;
      U[i_max][i_max] = cos(phi);
      U[j_max][j_max] = cos(phi);
      U[i_max][j_max] = -sin(phi);
      U[j_max][i_max] = sin(phi);
      U_T = U.transpose();
      A = U_T * A;
      A = A * U;
      X = X * U;
      k += 1;
   }

   Matrix Lambda = Matrix(1, n);
   for (int i = 0; i < n; ++i){
      Lambda[0][i] = A[i][i];
   }

   Lambda.print();
   cout << k << "\n";
}
