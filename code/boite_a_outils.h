#ifndef BOITE_A_OUTILS
#define BOITE_A_OUTILS

#include <vector>

#include "maillage.h"
#include "triangle.h"

using namespace std;

// On définit l'addition de 2 vecteur de double (ils doivent être de même
// taille).
vector<double> operator+(vector<double> A, vector<double> B);

// On définit la soustraction de 2 vecteur de double (ils doivent être de même
// taille).
vector<double> operator-(vector<double> A, vector<double> B);

// On définit le produit d'un vecteur de double avec un double.
vector<double> operator*(double a, vector<double> B);

// On définit le produit scalaire de 2 vecteur de double (ils doivent être de
// même taille).
double operator*(vector<double> A, vector<double> B);

// Question 3:
// Retourne f_eta (x, y).
double f_second_membre(double (*ug)(double), double (*ud)(double),
                       double (*ugpp)(double), double (*udpp)(double), double x,
                       double y, double a);

// Question 43.d:
// Calcule la prolongation du vecteur des noeuds intérieurs sur tous les
// noeuds.
vector<double> extend_vec(vector<double> V_int, int M, int N);

// Question 43.e:
// Calcule la restriction du vecteur des noeuds globaux sur tous les noeuds
// intérieurs.
vector<double> int_vec(vector<double> V_glb, int M, int N);

// Question 46.i:
// Retourne la norme L2 de la fonction v associée au vecteur V.
double norme_L2(vector<double> V, Maillage maille);

// Question 46.j:
// Retourne la norme L2 grad de la fonction v associée au vecteur V.
double norme_L2_grad(vector<double> V, Maillage maille);

// Question 44:
// Retourne le produit vectoriel A_eta * V.
vector<double> mat_vec(vector<double> V, Maillage maille);

// Question 45:
// Retourne le second membre B_eta du système linéaire A_eta * X = B_eta.
vector<double> scd_membre(double (*rhfs)(double, double), Maillage maille);

// Partie 4:
// Retourne une solution approchée du système linéaire A_eta * X = B_eta.
vector<double> inv_syst(vector<double> B_eta, Maillage maille,
                        int max_iteration);

#endif