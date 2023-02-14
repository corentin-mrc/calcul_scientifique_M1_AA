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
  int N = 10, M = N;
  double a = 1, b = a;
  Maillage maille(N, M, a, b);
  vector<double> B_eta = scd_membre(f_eta, maille);
  cout << "Le second membre est:" << endl;
  for (double i : B_eta)
    cout << i << endl;
  vector<double> w_eta_h = inv_syst(B_eta, maille, 10);
  cout << endl << "La solution approchée est:" << endl;
  for (double i : w_eta_h)
    cout << i << endl;
  cout << endl << "Les erreurs sont:" << endl;
  vector<double> erreur = erreurs(u_eta, w_eta_h, maille);
  for (double err : erreur)
    cout << err << endl;
  return 0;
}