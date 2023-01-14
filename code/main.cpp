#include <iostream>
#include <vector>

#include "constantes.h"
#include "noeud.h"

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
// Retourne une matrice K * 3 dont la l-ème ligne contient les trois numéros
// globaux des sommets du triangle T_l.
vector<vector<int>> maillageTR(int N, int M) {
  // Génère la matrice des triangles.
  vector<vector<int>> triangles;
  for (int j = 0; j < M; j++) {
    vector<int> triangle1;
    vector<int> triangle2;
    for (int i = 0; i < N; i++) {
      // Il y a 2 configurations possibles en fonction de la position du
      // rectangle du maillage qui nous intéresse. Dans chaque cas, il
      // faut déterminer les sommets du rectangle et les placer dans un
      // certain ordre.
      int a = Noeud(i, j).numgb(N, M);
      int b = Noeud(i + 1, j).numgb(N, M);
      int c = Noeud(i, j + 1).numgb(N, M);
      int d = Noeud(i + 1, j + 1).numgb(N, M);
      if (i % 2 == j % 2) {
        triangle1 = {a, c, d};
        triangle2 = {a, b, d};
      } else {
        triangle1 = {a, b, c};
        triangle2 = {b, c, d};
      }
      triangles.push_back(triangle1);
      triangles.push_back(triangle2);
    }
  }
  return triangles;
}

int main(void) {
  // Un petit test pour les noeuds.
  int N = 5;
  int M = 3;
  Noeud n(1, 2, 3.4, 3.4);
  cout << n.numgb(N, M) << endl;
  cout << n.numint(N, M) << endl;
  cout << n.num_gb_int(N, M, 21) << endl;
  cout << n.num_int_gb(N, M, 6) << endl;

  // Un petit test pour la triangulation.
  vector<vector<int>> tri = maillageTR(N, M);
  cout << "Test de la triangulation pour N = " << N << " et M = " << M << endl;
  for (vector<int> triangle : tri) {
    for (int i : triangle) {
      cout << i << " ";
    }
    cout << endl;
  }
  cout << "Il y a " << tri.size() << " triangles." << endl;

  return 0;
}
