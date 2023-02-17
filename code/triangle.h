#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>

#include "donnees_du_probleme.h"
#include "noeud.h"

using namespace std;

class Triangle {
  // sommets du triangle.
  vector<Noeud> noeuds;

public:
  // Constructeur par défaut.
  Triangle(void);
  // Constructeurs.
  Triangle(Noeud n0, Noeud n1, Noeud n2);
  Triangle(vector<Noeud> noeuds);

  // Getter.
  vector<Noeud> get_noeuds(void);

  // Question 31:
  // Retourne la matrice B_T du triangle.
  vector<vector<double>> calc_mat_BT(void);

  // Retourne la matrice B_T du triangle avec une permutation des sommets.
  vector<vector<double>> calc_mat_BT(vector<int> permut);

  // Retourne l'inverse de la matrice B_T.
  vector<vector<double>> inv_mat_BT(void);

  // Question 37:
  vector<vector<double>> diff_terme(void);
  vector<vector<double>> convect_terme(void);
  vector<vector<double>> react_terme(void);
};

#endif