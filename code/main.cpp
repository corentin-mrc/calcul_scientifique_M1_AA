#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "boite_a_outils.h"
#include "donnees_du_probleme.h"
#include "maillage.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

double epsilon = 1;
double gama = 1;
double lambda = 1;
double a = 1;
double b = 1;
int N = 12;
int M = 12;
double det = abs(2 * (a / N) * (b / M));

double u_g(double y) { return sin(M_PI * y); }
double u_gpp(double y) { return -M_PI * M_PI * sin(M_PI * y); }

double u_d(double y) { return 0; }
double u_dpp(double y) { return 0; }

double f_eta(double x, double y) {
  return f_second_membre(u_g, u_d, u_gpp, u_dpp, x, y);
}

int main(void) {
  int I = (N - 1) * (M - 1);
  Maillage maille;
  vector<double> B_eta = scd_membre(f_eta, maille);
  cout << "Le second membre est:" << endl;
  for (int k = 0; k < I && k < 30; k++)
    cout << B_eta[k] << endl;
  // w_eta_h est la solution approchée.
  vector<double> w_eta_h = inv_syst(B_eta, maille, 5);
  cout << endl << "La solution approchée est:" << endl;
  for (int k = 0; k < I && k < 30; k++)
    cout << w_eta_h[k] << endl;
  cout << endl << "Les erreurs sont:" << endl;
  vector<double> erreur = erreurs(u_eta, w_eta_h, maille);
  for (double err : erreur)
    cout << err << endl;
  vector<double> solution_exacte;
  for (int k = 0; k < I; k++) {
    vector<double> xy = maille.int_coord(k);
    solution_exacte.push_back(u_eta(xy[0], xy[1]));
  }
  vector<double> ecart = solution_exacte - w_eta_h;
  cout << endl << "sol exacte\tsol approchee\tecart" << endl;
  for (int k = 0; k < I && k < 30; k++)
    cout << solution_exacte[k] << "   \t" << w_eta_h[k] << "   \t" << ecart[k]
         << endl;
  cout << endl;

  // Ecriture des fichiers:
  ofstream file;

  // Ecriture des solutions exactes:
  file.open("solution_exacte.txt");
  for (double d : solution_exacte)
    file << d << endl;
  file.close();

  // Ecriture des solutions approchées:
  file.open("solution_approchee.txt");
  for (double d : w_eta_h)
    file << d << endl;
  file.close();
  
  
  vector<double> X0(B_eta.size(), 1);
  vector<double> AX0 = mat_vec(X0, maille);
  vector<double> AB_eta = mat_vec(B_eta, maille);
  cout << "A * X0  \tA * B_eta" << endl;
  for (int k = 0; k < I && k < 30; k++)
	  cout << AX0[k] << "   \t" << AB_eta[k] << endl;
  cout << endl;

  return 0;
}