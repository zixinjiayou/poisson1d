#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;



void poisson1d(ArrayXd& u, ArrayXd f,double alpha_l, double beta_l, double gamma_l,double alpha_h, double beta_h, double gamma_h, double h)
{
	// Find N. Declare local arrays.
	const int N = u.rows() ;     //taille de la premiere colonne de u

	MatrixXd A= Eigen::MatrixXd::Identity(N, N) ;
	MatrixXd A_inv;

		
	A = - A * 2;
	for (int i = 1; i <= N-2; i++) {
		A(i, i - 1) = 1;
		A(i, i + 1) = 1;
	}

	A_inv = A.inverse();

	VectorXd w = f * h * h;

	w(0) -= gamma_l * h / (alpha_l * h - beta_l);
	w(N-1) -= gamma_h * h / (alpha_h * h + beta_h);

	u = A_inv * w;
	// Calculate i=0 and i=N+1 values
	u(0) = (gamma_l * h - beta_l * u(1)) /
		(alpha_l * h - beta_l);
	u(N-1) = (gamma_h * h + beta_h * u(N-1)) /
		(alpha_h * h + beta_h);
};