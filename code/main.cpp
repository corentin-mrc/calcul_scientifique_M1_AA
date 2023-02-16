#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

#include "boite_a_outils.h"
#include "maillage.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

double u_g(double y) { return sin(M_PI * y); }
double u_gpp(double y) { return -M_PI * M_PI * sin(M_PI * y); }

double u_d(double y) { return 0; }
double u_dpp(double y) { return 0; }

double f_eta(double x, double y) {
  return f_second_membre(u_g, u_d, u_gpp, u_dpp, x, y, 1);
}

int main(void) {
  int N = 25, M = N;
  double a = 1, b = a;
  int I = (N - 1) * (M - 1);
  Maillage maille(N, M, a, b);
  vector<double> B_eta = scd_membre(f_eta, maille);
  cout << "Le second membre est:" << endl;
  for (int k = 0; k < I && k < 30; k++)
    cout << B_eta[k] << endl;
  // w_eta_h est la solution approchée.
  vector<double> w_eta_h = inv_syst(B_eta, maille, 25);
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
	  cout << solution_exacte[k] << "   \t" << w_eta_h[k] << "   \t" << ecart[k] << endl;
  cout << endl;
  return 0;
}