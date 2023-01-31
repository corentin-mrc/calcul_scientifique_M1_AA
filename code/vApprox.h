#ifndef VAPPROX_H
#define VAPPROX_H

#include <vector>

#include "maillage.h"
#include "noeud.h"

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
  // retourne la prolongation du vecteur des noeuds intérieurs sur tous les
  // noeuds
  void extendVec(void);

  // Question 43.e:
  // retourne la restriction du vecteur des noeuds globaux sur tous les noeuds
  // intérieurs
  void IntVec(void);
};

#endif