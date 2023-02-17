#ifndef BOITE_A_OUTILS
#define BOITE_A_OUTILS

#include <vector>

#include "donnees_du_probleme.h"
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
vector<double> operator*(double scalaire, vector<double> B);

// On définit le produit scalaire de 2 vecteur de double (ils doivent être de
// même taille).
double operator*(vector<double> A, vector<double> B);

// Retourne la plus grande valeur absolue du vecteur.
double max(vector<double> A);

// Question 3:
// Retourne f_eta (x, y).
double f_second_membre(double (*u_g)(double), double (*u_d)(double),
                       double (*u_gpp)(double), double (*u_dpp)(double),
                       double x, double y);

// Question 43.d:
// Calcule la prolongation du vecteur des noeuds intérieurs sur tous les
// noeuds.
vector<double> extend_vec(vector<double> V_int);

// Question 43.e:
// Calcule la restriction du vecteur des noeuds globaux sur tous les noeuds
// intérieurs.
vector<double> int_vec(vector<double> V_glb);

// Question 46.i:
// Retourne la norme L2 de la fonction v associée au vecteur V de taille I.
double norme_L2(vector<double> V, Maillage maille);

// Question 46.j:
// Retourne la norme L2 grad de la fonction v associée au vecteur V de taille I.
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

// Question 49:
// Retournes les trois erreurs relatives.
vector<double> erreurs(double (*sol_exa)(double, double),
                       vector<double> sol_appr, Maillage maille);

// TEMPORAIRE
double u_eta(double x, double y);

#endif