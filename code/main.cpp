#include <iostream>
#include <vector>

#include "boite_a_outils.h"
#include "maillage.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

int main(void) {
  // Un petit test pour la triangulation.
  int N = 5;
  int M = 4;
  double a = 10;
  double b = 10;
  Maillage maille(N, M, a, b);
  for (Triangle triangle : maille.get_triangulation()) {
    for (Noeud noeud : triangle.get_noeuds()) {
      cout << maille.num_gb_noeud(noeud) << " ";
    }
    cout << endl;
  }

  // Un petit test pour les surcharges.
  vector<double> victor = {1, 2, 3, 4, 5};
  vector<double> victoire = {6, 1, 0, -2, -12};
  vector<double> adrien = victor + victoire;
  double guy = victor * victoire;
  vector<double> pascal = 3.14 * victoire;
  for (int i = 0; i < 5; i++)
    cout << victor[i] << "\t+\t" << victoire[i] << "\t=\t" << adrien[i] << endl;
  cout << "produit scalaire : " << guy << endl;
  for (int i = 0; i < 5; i++)
    cout << "3.14\t*\t" << victoire[i] << "\t=\t" << pascal[i] << endl;
  cout << endl;

  // Un petit test pour les normes.
  vector<double> v = {1, 3, 2, 4, 5, 9, 0, 8, 6, 7, 4, 5};
  cout << norme_L2(v, maille) << endl;
  cout << norme_L2_grad(v, maille) << endl;

  // Un petit test pour le produit matrice vecteur.
  vector<double> v_test = mat_vec(v, maille);
  for (double i : v_test)
    cout << i << " ";
  cout << endl;

  return 0;
}