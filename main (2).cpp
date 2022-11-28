#include <iostream>
#include <Eigen/Dense>
#include <fstream> // pour pouvoir �crire dans des fichiers
#include <cmath>

using namespace Eigen;
using namespace std;


void poisson1d(ArrayXd& u, ArrayXd f,
	double alpha_l, double beta_l, double gamma_l,
	double alpha_h, double beta_h, double gamma_h, double dx);


int main()
{
	// d�claration de la fonction poisson1d.
	// voir poisson1d.cpp pour l�impl�mentation

	const int n = 100; // nombre de points de discr�tisation

	ArrayXd i = ArrayXd::LinSpaced(n + 1, 0, n);// indice varaint de 0 a n+1

	double l = 1.0; // intervalle [0,l]
	double h = l / (n + 1); // pas du maillage

	ArrayXd ucal (n + 1); // solution 
	ArrayXd solexacte(n + 1); // solution exacte
	ArrayXd sndmbre(n + 1); // second membre


	//exemple de solution exacte et second membre
	double coeff = 250.;

	sndmbre = -(2 * coeff * (1. + coeff * ((h * i - 0.5) * (i * h - 0.5))) -
		8 * coeff * coeff * (i * h - 0.5) * (i * h - 0.5)) /
		pow((1. + coeff * (i * h - 0.5) * (i * h - 0.5)), 3);

	solexacte = 1. / (1. + coeff * (i * h - 0.5 * (i * h - 0.5)));

	// constantes pour les conditions aux limites
	double alpha_l = 1;
	double beta_l = 0;
	double gamma_l = 1. / (1. + coeff * 0.25);
	double alpha_h = 1;
	double beta_h = 0;
	double gamma_h = gamma_l;


	// appel de la fonction poisson1d pour r�soudre le probl�me


	poisson1d(ucal, sndmbre,alpha_l, beta_l, gamma_l,alpha_h, beta_h, gamma_h, h);


	// calcul de l�erreur entre solution exacte et solution calcul�e

	cout << "root-mean-square error:" << endl
		<< "2nd order approximation: " << sqrt(sqrt(abs(solexacte - ucal)).mean())
			<< endl;

	ofstream sortiefichier("lapl1d.dat");

	cout << endl;
	cout << "creation du fichier �lapl1d.dat� " << endl;
	cout << endl;

	sortiefichier << i * h << " " << ucal << " "<< solexacte << endl;
	sortiefichier.close(); // fermeture du fichier
	cout << "***" << endl;
	cout << "pour visualiser la solution exacte et la solution que vous "
		<< "venez de calculer, utilisez la ligne de commande suivante "
		<< " sous �gnuplot� : " << endl;
	cout << "plot �lapl1d.dat� using 1:2 title �ucal�,"
		<< "�lapl1d.dat� using 1:3 with lines title �uex�" << endl;
	cout << "***" << endl;
	return 0;
}



