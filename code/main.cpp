#include <iostream>
#include <vector>

#include "constantes.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

// Question 3:
// Retourne f_eta (x, y).
double fsecond_membre(double (*ug)(double), double (*ud)(double),
                      double (*ugpp)(double), double (*udpp)(double), double x,
                      double y) {
  return (-epsilon * (a - x) * ugpp(y) - epsilon * (a + x) * udpp(y) -
          gama * ug(y) + gama * ud(y) + lambda * (a - x) * ug(y) +
          lambda * (a + x) * ud(y)) /
         (2 * a);
}

// Question 9:
// Retourne une subdivision uniforme de [-a, a] en N + 1 points.
vector<double> Subdiv(int N) {
  vector<double> xi;
  for (int i = 0; i <= N; i++)
    xi.push_back(-a + (2. * i * a) / N);
  return xi;
}

// Question 20.
// Retourne un tableau de triangle dont la l-ème ligne contient le triangle T_l.
vector<Triangle> maillageTR(int N, int M) {
  // Génère la matrice des triangles.
  vector<Triangle> triangles;
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < N; i++) {
      // Il y a 2 configurations possibles en fonction de la position du
      // rectangle du maillage qui nous intéresse. Dans chaque cas, il
      // faut déterminer les sommets du rectangle et les placer dans un
      // certain ordre.
      Noeud a(i, j);
      Noeud b(i + 1, j);
      Noeud c(i, j + 1);
      Noeud d(i + 1, j + 1);
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
  int N = 10;
  int M = 20;
  Noeud n(1, 2, 3.4, 3.4);
  cout << n.numgb(N, M) << endl;
  cout << n.numint(N, M) << endl;
  cout << n.num_gb_int(N, M, 21) << endl;
  cout << n.num_int_gb(N, M, 6) << endl;

  cout << endl;

  // Un petit test pour la triangulation.
  vector<Triangle> tri = maillageTR(N, M);
  cout << "Test de la triangulation pour N = " << N << " et M = " << M << endl;
  for (Triangle triangle : tri)
    triangle.affiche_sommets_glb(N, M);
  cout << "Il y a " << tri.size() << " triangles." << endl;

  return 0;
}
