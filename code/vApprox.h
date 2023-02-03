#ifndef VAPPROX_H
#define VAPPROX_H

#include <vector>

#include "maillage.h"
#include "noeud.h"
#include "triangle.h"

using namespace std;

class VApprox {
  vector<double> vInt;
  vector<double> vGlb;
  Maillage maille;

public:
  // Constructeur par défaut:
  VApprox(void);
  // Constructeur:
  VApprox(vector<double> vInt, Maillage maille);

  // Getters:
  vector<double> getvInt(void);
  vector<double> getvGlb(void);

  // Setter:
  void setvGlb(vector<double> vGlb);

  // Question 43.d:
  // Calcule la prolongation du vecteur des noeuds intérieurs sur tous les
  // noeuds.
  void extendVec(void);

  // Question 43.e:
  // Calcule la restriction du vecteur des noeuds globaux sur tous les noeuds
  // intérieurs.
  void IntVec(void);

  // Retourne la restriction du vecteur des noeuds globaux sur tous les noeuds
  // intérieurs.
  static vector<double> IntVec(vector<double> vGlb, Maillage maille);
  
  // Question 46.i:
  // Retourne la norme L2 de la fonction v associée au vecteur V
  double normL2(vector<Triangle> triangles);
  
  // Question 46.j:
  // Retourne la norme L2 grad de la fonction v associée au vecteur V
  double normL2Grad(vector<Triangle> triangles);
};

#endif