#include <iostream>
#include <vector>

#include "constantes.h"
#include "maillage.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

// Question 3:
// Retourne f_eta (x, y).
double fsecond_membre(double (*ug)(double), double (*ud)(double),
                      double (*ugpp)(double), double (*udpp)(double), double x,
                      double y, double a) {
  return (epsilon * (a - x) * ugpp(y) + epsilon * (a + x) * udpp(y) +
          gama * ug(y) - gama * ud(y) - lambda * (a - x) * ug(y) -
          lambda * (a + x) * ud(y)) /
         (2 * a);
}

// Question 9:
// Retourne une subdivision uniforme de [-a, a] en N + 1 points.
vector<double> Subdiv(double a, int N) {
  vector<double> xi;
  for (int i = 0; i <= N; i++)
    xi.push_back(-a + (2 * i * a) / N);
  return xi;
}

// Question 20.
// Retourne un tableau de triangle dont la l-ème ligne contient le triangle T_l.
vector<Triangle> maillageTR(Maillage maille) {
  // Génère la matrice des triangles.
  vector<Triangle> triangles;
  for (int j = 0; j < maille.getM(); j++) {
    for (int i = 0; i < maille.getN(); i++) {
      // Il y a 2 configurations possibles en fonction de la position du
      // rectangle du maillage qui nous intéresse. Dans chaque cas, il
      // faut déterminer les sommets du rectangle et les placer dans un
      // certain ordre.
      Noeud a(i, j, maille);
      Noeud b(i + 1, j, maille);
      Noeud c(i, j + 1, maille);
      Noeud d(i + 1, j + 1, maille);
      if (i % 2 == j % 2) {
        triangles.push_back(Triangle(a, c, d));
        triangles.push_back(Triangle(a, b, d));
      } else {
        triangles.push_back(Triangle(a, b, c));
        triangles.push_back(Triangle(b, c, d));
      }
    }
  }
  return triangles;
}

int main(void) {
  // Un petit test pour les noeuds.
  int N = 5;
  int M = 4;
  Maillage m(N, M, 10, 10);
  Noeud n(1, 2, m);
  cout << n.numgb(m) << endl;
  cout << n.numint(m) << endl;
  cout << n.num_gb_int(m, 21) << endl;
  cout << n.num_int_gb(m, 6) << endl;

  cout << endl;

  // Un petit test pour la triangulation.
  vector<Triangle> triangulation = maillageTR(m);
  cout << "Test de la triangulation pour N = " << N << " et M = " << M << endl;
  for (Triangle triangle : triangulation)
    triangle.affiche_sommets_glb(m);
  cout << "Il y a " << triangulation.size() << " triangles." << endl;

  cout << endl;
  Triangle t(Noeud(0, 0, m), Noeud(1, 0, m), Noeud(0, 1, m));
  cout << t.DetMatBT() << endl;

  return 0;
}