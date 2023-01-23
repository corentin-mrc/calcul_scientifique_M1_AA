#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>

#include "maillage.h"
#include "noeud.h"

using namespace std;

class Triangle {
  Noeud n0;
  Noeud n1;
  Noeud n2;

public:
  // Constructeur par défaut:
  Triangle(void);
  // Constructeur:
  Triangle(Noeud n0, Noeud n1, Noeud n2);

  // Question 31:
  // Retourne la matrice B_T du triangle.
  vector<vector<double>> CalcMatBT(void);

  // Retourne le déterminant de la matrice B_T:
  double DetMatBT(void);

  // Retourne l'inverse de la matrice B_T:
  vector<vector<double>> InvMatBT(void);

  // Question 37:
  vector<vector<double>> DiffTerm(void);

  vector<vector<double>> ConvectTerm(void);

  vector<vector<double>> ReacTerm(void);
  
  // retourne les noeuds
  vector<Noeud> noeuds(void);

  // Affiche les numéros globaux des sommets du triangle.
  void affiche_sommets_glb(Maillage maille);
};

#endif