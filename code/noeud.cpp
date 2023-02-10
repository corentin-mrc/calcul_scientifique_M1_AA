#include "noeud.h"

Noeud::Noeud(void) {}

Noeud::Noeud(double x, double y) {
  this->x = x;
  this->y = y;
}

double Noeud::get_x(void) { return x; }

double Noeud::get_y(void) { return y; }