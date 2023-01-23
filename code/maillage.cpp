#include <vector>

#include "maillage.h"

Maillage::Maillage(void) {}

Maillage::Maillage(int N, int M, double a, double b) {
  this->N = N;
  this->M = M;
  this->a = a;
  this->b = b;
}

int Maillage::getN(void) { return N; }

int Maillage::getM(void) { return M; }

double Maillage::geta(void) { return a; }

double Maillage::getb(void) { return b; }

vector<vector<double>> Maillage::Subdiv(void) {
  vector<double> xi;
  for (int i = 0; i <= N; i++)
    xi.push_back(-a + (2 * i * a) / N);
  vector<double> yi;
  for (int j = 0; j <= M; j++)
    yi.push_back(-b + (2 * j * b) / M);
  return {xi, yi};
}
