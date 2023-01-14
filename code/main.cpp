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

int main(void) {
  // Un petit test pour les noeuds.
  int N = 5;
  int M = 4;
  Noeud n(1, 2, 3.4, 3.4);
  cout << n.numgb(N, M) << endl;
  cout << n.numint(N, M) << endl;
  cout << n.num_gb_int(N, M, 21) << endl;
  cout << n.num_int_gb(N, M, 6) << endl;
  return 0;
}
