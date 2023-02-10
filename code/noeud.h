#ifndef NOEUD_H
#define NOEUD_H

class Noeud {
  // Coordonnées du noeud
  double x, y;

public:
  // Constructeur par défaut.
  Noeud(void);
  // Constructeur.
  Noeud(double x, double y);

  // Getters.
  double get_x(void);
  double get_y(void);
};

#endif