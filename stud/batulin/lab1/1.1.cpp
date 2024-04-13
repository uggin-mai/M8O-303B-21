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

   void inverse(Matrix &L, Matrix &U, double Z[])
   {
      Matrix E = Matrix(rows, cols, 1);
      for (int k = 0; k < rows; ++k){
         for (int i = 0; i < cols; ++i){
            double s = 0;
            for (int j = 0; j < i; ++j){
               s += L[i][j] * Z[j];
            }
            Z[i] = E[i][k] - s;
         }

         for (int i = rows - 1; i > -1; --i){
            double s = 0;
            for (int j = i + 1; j < rows; ++j){
               s += U[i][j] * a[j][k];
            }
            a[i][k] = (Z[i] - s) / U[i][i];
         }
      }
   }

   double det()
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
            d += pow(-1, k + 2) * a[0][k] * M.det();
         }
      return d;
   }

	double* operator[](int i)
	{
      return a[i];
   	}

   Matrix& operator*=(Matrix& m)
   {    
      if (cols != m.rows)
      {
      throw "Wait. That's illegal.";
      }
      Matrix temp(m.rows, m.cols);
      for (int i = 0; i < temp.rows; ++i) 
      {
         for (int j = 0; j < temp.cols; ++j) 
         {
               for (int k = 0; k < cols; ++k) 
               {
                  temp.a[i][j] += (a[i][k] * m.a[k][j]);
               }
         }
      }
      return (*this = temp);
   }
};

void LU(Matrix &A, Matrix &B, Matrix &L, Matrix &U, Matrix M[])
{
   int n = A.rows;
   L = Matrix(n, n, 1);

   for(int i = 0;i < n; i++)
    for(int j = 0;j < n; j++)
        U[i][j]=A[i][j];

   for (int k = 0; k < n - 1; k++){
      M[k] = Matrix(n, n, 1);
      for (int i = k + 1; i < n; i++){
         if (U[k][k] == 0){
            int j = k + 1;
            while (U[j][j] == 0 and j < n){
               j += 1;
            }
            if (j == n){
               break;
            }
            U.swapRows(k, j);
            B.swapRows(k, j);
         }
         M[k][i][k] = U[i][k] / U[k][k];
         for (int j = k; j < n; ++j){
            U[i][j] -= M[k][i][k] * U[k][j];
         }
      }
      L *=M[k];
   }
}

/*
1 2 -2 6 -3 -5 14 13 1 2 -2 -2 -2 -4 5 10
24 41 0 20
*/

int main()
{
   int n = 4;
   Matrix A = Matrix(n, n);
   Matrix U = Matrix(n, n);
   Matrix B = Matrix(n, 1);
   A.scan();
   B.scan();

   Matrix M[n - 1];
   Matrix L = Matrix(n, n, 1);

   LU(A,B,L,U,M);

   double Z[n];
   for (int i = 0; i < n; ++i){
      double s = 0;
      for (int j = 0; j < i; ++j){
         double tmp = L[i][j] * Z[j];
         s += tmp;
      }
      Z[i] = B[i][0] - s;
   }

   double X[n];
   for (int i = n - 1; i > -1; --i){
      double s = 0;
      for (int j = i + 1; j < n; ++j){
         s += U[i][j] * X[j];
      }
      X[i] = (Z[i] - s) / U[i][i];
   }

   for (int i = 0; i < n; ++i){
      cout << X[i] << " ";
   }
   cout << "\n\n";

   cout << U.det() << "\n\n";

   Matrix Ai = Matrix(n, n);
   Ai.inverse(L,U,Z);
   Ai.print();
}
