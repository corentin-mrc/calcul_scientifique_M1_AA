#include <vector>

#include "noeud.h"

using namespace std;

Noeud::Noeud(void) {}

Noeud::Noeud(int i, int j) {
  this->i = i;
  this->j = j;
}

Noeud::Noeud(int i, int j, Maillage maille) {
  this->i = i;
  this->j = j;
  this->x = 2 * i * maille.geta() / maille.getN();
  this->y = 2 * j * maille.getb() / maille.getM();
}

Noeud::Noeud(int i, int j, double x, double y) {
  this->i = i;
  this->j = j;
  this->x = x;
  this->y = y;
}

double Noeud::getx(void) { return x; }

double Noeud::gety(void) { return y; }

int Noeud::numgb(Maillage maille) { return (maille.getN() + 1) * j + i; }

void Noeud::invnumgb(Maillage maille, int s) {
  i = s % (maille.getN() + 1);
  j = s / (maille.getN() + 1);
}

int Noeud::numint(Maillage maille) {
  return (maille.getN() - 1) * (j - 1) + (i - 1);
}

void Noeud::invumint(Maillage maille, int k) {
  i = k % (maille.getN() - 1) + 1;
  j = k / (maille.getN() - 1) + 1;
}

int Noeud::num_int_gb(Maillage maille, int k) {
  invumint(maille, k);
  return numgb(maille);
}

int Noeud::num_gb_int(Maillage maille, int s) {
  invnumgb(maille, s);
  return numint(maille);
}