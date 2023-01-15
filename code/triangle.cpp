#include <iostream>

#include "triangle.h"

using namespace std;

Triangle::Triangle(void) {}

Triangle::Triangle(Noeud n1, Noeud n2, Noeud n3) {
  this->n1 = n1;
  this->n2 = n2;
  this->n3 = n3;
}

void Triangle::affiche_sommets_glb(int N, int M) {
  cout << n1.numgb(N, M) << " " << n2.numgb(N, M) << " " << n3.numgb(N, M)
       << endl;
}