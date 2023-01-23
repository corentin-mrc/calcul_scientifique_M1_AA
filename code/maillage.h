#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <vector>

using namespace std;

class Maillage {
  int N;
  int M;

  double a;
  double b;

public:
  // Constructeur par défault.
  Maillage(void);
  // Constructeurs.
  Maillage(int N, int M, double a, double b);

  // Getters.
  int getN(void);
  int getM(void);
  double geta(void);
  double getb(void);

  // Question 9:
  // Retourne une subdivision uniforme de [-a, a] en N + 1 points et de [-b, b]
  // en M + 1 points.
  vector<vector<double>> Subdiv(void);
};

#endif