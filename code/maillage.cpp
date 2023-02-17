#include <cmath>

#include "maillage.h"

Maillage::Maillage(void) { this->init_maillage_TR(); }

vector<Triangle> Maillage::get_triangulation(void) { return triangulation; }

vector<double> Maillage::sub_div(double largeur, int nb_divisions) {
  vector<double> xi;
  for (int i = 0; i <= nb_divisions; i++)
    xi.push_back(-largeur + (2 * i * largeur) / nb_divisions);
  return xi;
}

int Maillage::num_gb(int i, int j) { return (N + 1) * j + i; }

vector<int> Maillage::inv_num_gb(int s) {
  int i = s % (N + 1);
  int j = s / (N + 1);
  return {i, j};
}

int Maillage::num_int(int i, int j) { return (N - 1) * (j - 1) + (i - 1); }

vector<int> Maillage::inv_num_int(int k) {
  int i = k % (N - 1) + 1;
  int j = k / (N - 1) + 1;
  return {i, j};
}

int Maillage::num_int_gb(int k) {
  vector<int> ij = this->inv_num_int(k);
  return this->num_gb(ij[0], ij[1]);
}

int Maillage::num_gb_int(int s) {
  vector<int> ij = this->inv_num_gb(s);
  return this->num_int(ij[0], ij[1]);
}

bool Maillage::est_sur_le_bord(Noeud noeud) {
  int i = round(N * (noeud.get_x() + a) / (2 * a));
  int j = round(M * (noeud.get_y() + b) / (2 * b));
  return i == 0 || j == 0 || i == N || j == M;
}

int Maillage::num_gb_noeud(Noeud noeud) {
  int i = round(N * (noeud.get_x() + a) / (2 * a));
  int j = round(M * (noeud.get_y() + b) / (2 * b));
  return this->num_gb(i, j);
}

int Maillage::num_int_noeud(Noeud noeud) {
  int i = round(N * (noeud.get_x() + a) / (2 * a));
  int j = round(M * (noeud.get_y() + b) / (2 * b));
  return this->num_int(i, j);
}

vector<double> Maillage::int_coord(int k) {
  vector<int> ij = this->inv_num_int(k);
  return {(2 * ij[0] - N) * a / N, (2 * ij[1] - M) * b / M};
}

void Maillage::init_maillage_TR(void) {
  // On génère la matrice des noeuds.
  vector<vector<Noeud>> noeuds;
  for (int i = 0; i <= N; i++) {
    vector<Noeud> colonne;
    double x = (2 * i - N) * a / N;
    for (int j = 0; j <= M; j++) {
      double y = (2 * j - M) * b / M;
      colonne.push_back(Noeud(x, y));
    }
    noeuds.push_back(colonne);
  }
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < N; i++) {
      Noeud SO = noeuds[i][j];
      Noeud SE = noeuds[i + 1][j];
      Noeud NO = noeuds[i][j + 1];
      Noeud NE = noeuds[i + 1][j + 1];
      // Il y a 2 configurations possibles en fonction de la position du
      // rectangle du maillage qui nous intéresse. Dans chaque cas, il
      // faut déterminer les sommets du rectangle et les placer dans un
      // certain ordre.
      if ((i ^ j) & 1) {
        triangulation.push_back(Triangle(SO, SE, NO));
        triangulation.push_back(Triangle(SE, NO, NE));
      } else {
        triangulation.push_back(Triangle(SO, NO, NE));
        triangulation.push_back(Triangle(SO, SE, NE));
      }
    }
  }
}
