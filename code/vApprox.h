#ifndef VAPPROX_H
#define VAPPROX_H


#include <vector>

#include "maillage.h"
#include "noeud.h"

using namespace std;

class VApprox {
  Maillage maille;
  vector<double> vInt;
  vector<double> vGlb;

public:
  // Constructeur par défaut:
  VApprox(void);
  // Constructeur:
  VApprox(Maillage maille, vector<double> vInt);

  // Getters
  vector<double> getvInt(void);
  vector<double> getvGlb(void);
  
  // Question 43.d:
  // retourne la prolongation du vecteur des noeuds intérieurs sur tous les noeuds
  void extendVec(void);
  
  // Question 43.e:
  // retourne la restriction du vecteur des noeuds globaux sur les noeuds intérieurs
  void IntVec(void);
};



#endif