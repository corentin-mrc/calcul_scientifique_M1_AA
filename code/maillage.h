#ifndef MAILLAGE_H
#define MAILLAGE_H

#include "noeud.h"
#include "triangle.h"
#include <vector>

using namespace std;

class Maillage {
  // Ensemble de travail.
  double a, b;
  // Subdivisions.
  int N, M;
  vector<Triangle> triangulation;

public:
  // Constructeur par défault.
  Maillage(void);
  // Constructeur.
  Maillage(int N, int M, double a, double b);

  // Getters.
  int get_N(void);
  int get_M(void);
  double get_a(void);
  double get_b(void);
  vector<Triangle> get_triangulation(void);

  // Question 9:
  // Retourne une subdivision uniforme de [-a, a] en N + 1 points et de [-b, b]
  // en M + 1 points.
  static vector<double> sub_div(int a, int N);

  // Question 11:
  // Retourne le numéro global associé aux indices i et j.
  int num_gb(int i, int j);

  // Question 13:
  // Retourne les indices i et j à partir du numéro global s.
  vector<int> inv_num_gb(int s);

  // Question 14:
  // Retourne le numéro intérieur associé aux indices i et j.
  int num_int(int i, int j);

  // Question 16:
  // Retourne les indices i et j à partir du numéro intérieur k.
  vector<int> inv_num_int(int k);

  // Question 17:
  // Retourne le numéro global à partir du numéro intérieur k.
  int num_int_gb(int k);

  // Question 18:
  // Retourne le numéro intérieur à partir du numéro global s.
  int num_gb_int(int s);

  // Retourne le numéro global d'un noeud dans le maillage.
  int num_gb_noeud(Noeud noeud);

  // Retourne le numéro intérieur d'un noeud dans le maillage.
  int num_int_noeud(Noeud noeud);

  // Vérifie si le noeud est sur le bord.
  bool est_sur_le_bord(Noeud noeud);

  // Donne les coordonnées x et y à partir du numéro intérieur du noeud.
  vector<double> int_coord(int k);

  // Question 20:
  // Initialise le tableau de triangle dont la l-ème ligne contient le triangle
  // T_l.
  void init_maillage_TR(void);
};

#endif