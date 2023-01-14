#include "noeud.h"

Noeud::Noeud(int i, int j) {
  this->i = i;
  this->j = j;
}

Noeud::Noeud(int i, int j, double x, double y) {
  this->i = i;
  this->j = j;
  this->x = x;
  this->y = y;
}

int Noeud::numgb(int N, int M) { return (N + 1) * j + i; }

void Noeud::invnumgb(int N, int M, int s) {
  i = s % (N + 1);
  j = s / (N + 1);
}

int Noeud::numint(int N, int M) { return (N - 1) * (j - 1) + (i - 1); }

void Noeud::invumint(int N, int M, int k) {
  i = k % (N - 1) + 1;
  j = k / (N - 1) + 1;
}

int Noeud::num_int_gb(int N, int M, int k) {
  invumint(N, M, k);
  return numgb(N, M);
}

int Noeud::num_gb_int(int N, int M, int s) {
  invnumgb(N, M, s);
  return numint(N, M);
}