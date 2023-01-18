#ifndef NOEUD_H
#define NOEUD_H

#include "maillage.h"

class Noeud {
  // Indice dans le maillage global.
  int i;
  int j;

  // Coordonnées.
  double x;
  double y;

public:
  // Constructeur par défaut:
  Noeud(void);
  // Constructeurs:
  Noeud(int i, int j);
  Noeud(int i, int j, Maillage maille);
  Noeud(int i, int j, double x, double y);

  // Getters
  double getx(void);
  double gety(void);

  // Question 11:
  // Retourne le numéro global d'un noeud.
  int numgb(Maillage maille);

  // Question 13:
  // Détermine les indices i et j à partir du numéro global du noeud.
  void invnumgb(Maillage maille, int s);

  // Question 14:
  // Retourne le numéro intérieur d'un noeud.
  int numint(Maillage maille);

  // Question 16:
  // Determine les indices i et j à partir du numéro intérieur du noeud.
  void invumint(Maillage maille, int k);

  // Question 17:
  // Retourne le numéro global d'un noeud à partir de son numéro intérieur
  // Met à jour les coordonnées en fonction du numéro intérieur fourni.
  int num_int_gb(Maillage maille, int k);

  // Question 18:
  // Retourne le numéro intérieur d'un noeud à partir de son numéro global.
  // Met à jour les coordonnées en fonction du numéro global fourni.
  int num_gb_int(Maillage maille, int s);
};

#endif