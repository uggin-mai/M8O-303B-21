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

	double* operator[](int i)
	{
      return a[i];
   	}
};


double determinant(Matrix M)
{
	double det = 1;
   int p = 1;
	for (int i=0; i < M.rows; i++)
   {
		det *= pow(-1, p) * M[i][i];
      p++;
	}
	return det;
}

/*
0 -11 -9 5 -15 -2 -8 11 -3 6 -15 4 3 6 0
-122 -48 -14 -50 42
*/

int main()
{

   int n = 5;
   Matrix A = Matrix(n, 3);
   double B[n] = {-122, -48, -14, -50, -42};
   A.scan();
   A.print();
   //B.scan();

   // прямой ход прогонки
   double P[n], Q[n], X[n];
   P[0] = -A[0][2] / A[0][1];
   Q[0] = B[0] / A[0][1];
   for (int i = 1; i < n; ++i){
      P[i] = -A[i][2] / (A[i][1] + A[i][0] * P[i - 1]);
      Q[i] = (B[i] - A[i][0] * Q[i - 1]) / (A[i][1] + A[i][0] * P[i - 1]);
   }

   // обратный ход прогонки
   X[-1] = Q[-1];
   for (int i = n - 2; i > -1; --i){
      X[i] = P[i] * X[i + 1] + Q[i];
   }
   for (int i = 0; i < n; ++i){
      cout << X[i] << " ";
   }
}

    for(int i = 0; i < n; ++i) {
        if(i == 0) {
            P[i] = -A[i][i+1] / A[i][i];
            Q[i] = d[i] / A[i][i];
        } else if(i == n - 1) {
            P[i] = 0;
            Q[i] = (d[i] - A[i][i - 1] * Q[i - 1]) / (A[i][i] + A[i][i - 1] * P[i - 1]);
        } else {
            P[i] = -A[i][i+1] / (A[i][i] + A[i][i - 1] * P[i - 1]);
            Q[i] = (d[i] - A[i][i - 1] * Q[i - 1]) / (A[i][i] + A[i][i - 1] * P[i - 1]);
        }
    }

    for(int i = n - 1; i >= 0; --i) {
        if(i == n - 1) x[i] = Q[i];
        else {
            x[i] = P[i] * x[i + 1] + Q[i];
        }
    }
