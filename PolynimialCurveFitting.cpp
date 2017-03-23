// PolynimialCurveFitting.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>

using namespace std;


int main()
{
	int i, j, k, n, N;
	cout.precision(4);
	cout.setf(ios::fixed);
	cout << "\nEnter the number of data pairs to be entered:\n";
	cin >> N;
	
	vector<double> x(N), y(N);
	cout << "\nEnter the x-axis values:\n";
	for (i = 0; i < N; ++i) {
		cin >> x[i];
	}
	cout << "\nEnter the y-axis values:\n";
	for (i = 0; i < N; ++i) {
		cin >> y[i];
	}
	cout << "\nEnter the degree of Polynomial you want to use for the fit:\n";
	cin >> n;
	/*if (n % 10 != 0) {
		cout << "\nError: Invalid input! Must give an integer.\n";
		cout << "\nEnter the number of data pairs to be entered:\n";
		cin >> n;
	}*/

	// X is used to store the values of sigma(xi^0) = N, sigma(xi), sigma(xi^2)...sigma(xi^2n)
	vector<double> X(2*n + 1); 

	for (i = 0; i < 2 * n + 1; ++i) {
		X[i] = 0;
		for (j = 0; j < N; ++j) {
			X[i] += pow(x[j], i);
		}
	}
	//B is the Normal matrix (argmented) that will sotre the equations
	vector<vector<double>> B;
	B.resize(n + 1, vector<double>(n + 2));
	//a is for the value of the final coefficients, which we are about to obtain
	vector<double> a(n+1);
	for (i = 0; i <= n; ++i) {
		for (j = 0; j <= n; ++j) {
			B[i][j] = X[i + j];
		}
	}

	//Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
	//The matrix of sigma(yi), sigma(xi*yi), sigma(xi^2*yi)...sigma(xi^n*yi)
	vector<double> Y(n + 1);
	for (i = 0; i < n + 1; ++i) {
		Y[i] = 0;
		for (j = 0; j < N; ++j) {
			Y[i] += pow(x[j], i)*y[j];
		}
	}
	for (i = 0; i <= n; ++i) {
		B[i][n + 1] = Y[i];
	}
	//for n degree we get n+1 equations (Gaussian Elimination)
	n = n + 1;
	cout << "\nThe Norml (Augmented Matrix) is as following: \n";
	//To print out the Normal Augmented matrix
	for (i = 0; i < n; ++i) {
		for (j = 0; j <= n; ++j) {
			cout << B[i][j] << setw(16);
		}
		cout << "\n";
	}
	for (i = 0; i < n; ++i) {
		for (k = i + 1; k < n; ++k) {
			if (B[i][i] < B[k][i]) {
				for (j = 0; j <= n; ++j) {
					double temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}
			}
		}
	}
	for (i = 0; i < n - 1; ++i) { //loop to form the Guass Elimination
		for (k = i + 1; k < n; ++k) {
			double t = B[k][i] / B[i][i];
			for (j = 0; j <= n; ++j) {
				B[k][j] = B[k][j] - t * B[i][j]; 
				// make the elements below the pivot elements equal to zero or eliminate the variables
			}
		}
	}
	for (i = n - 1; i >= 0; --i) {
		a[i] = B[i][n];
		for (j = 0; j < n; ++j) {
			if (j != i) {
				a[i] = a[i] - B[i][j] * a[j];
			}
		}
		a[i] = a[i] / B[i][i];
	}
	cout << "\nThe values of the coefficients are as follows:\n";
	for (i = 0; i < n; ++i) {
		cout << "x^" << i << " = " << a[i] << endl;
	}
	cout << "\nHence the fitted Polynomila is given by: \n";
	cout << "y = ";
	for (i = 0; i < n; ++i) {
		cout << " + (" << a[i] << ")x^" << i;
	}
	//cout << "\n";
	return 0;
}

