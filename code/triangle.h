#include "noeud.h"

class Triangle {
  Noeud n1;
  Noeud n2;
  Noeud n3;

public:
  // Constructeur par défaut:
  Triangle(void);
  // Constructeur:
  Triangle(Noeud n1, Noeud n2, Noeud n3);

  // Affiche les numéros globaux des sommets du triangle.
  void affiche_sommets_glb(int N, int M);
};