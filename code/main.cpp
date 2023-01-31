#include <iostream>
#include <vector>

#include "constantes.h"
#include "maillage.h"
#include "noeud.h"
#include "triangle.h"
#include "vApprox.h"

using namespace std;

// Question 3:
// Retourne f_eta (x, y).
double fsecond_membre(double (*ug)(double), double (*ud)(double),
                      double (*ugpp)(double), double (*udpp)(double), double x,
                      double y, double a) {
  return (epsilon * (a - x) * ugpp(y) + epsilon * (a + x) * udpp(y) +
          gama * ug(y) - gama * ud(y) - lambda * (a - x) * ug(y) -
          lambda * (a + x) * ud(y)) /
         (2 * a);
}

// Question 9:
// Retourne une subdivision uniforme de [-a, a] en N + 1 points.
vector<double> Subdiv(double a, int N) {
  vector<double> xi;
  for (int i = 0; i <= N; i++)
    xi.push_back(-a + (2 * i * a) / N);
  return xi;
}

// Question 20:
// Retourne un tableau de triangle dont la l-ème ligne contient le triangle T_l.
vector<Triangle> maillageTR(Maillage maille) {
  // Génère la matrice des triangles.
  vector<Triangle> triangles;
  for (int j = 0; j < maille.getM(); j++) {
    for (int i = 0; i < maille.getN(); i++) {
      // Il y a 2 configurations possibles en fonction de la position du
      // rectangle du maillage qui nous intéresse. Dans chaque cas, il
      // faut déterminer les sommets du rectangle et les placer dans un
      // certain ordre.
      Noeud a(i, j, maille);
      Noeud b(i + 1, j, maille);
      Noeud c(i, j + 1, maille);
      Noeud d(i + 1, j + 1, maille);
      if (i % 2 == j % 2) {
        triangles.push_back(Triangle(a, c, d));
        triangles.push_back(Triangle(a, b, d));
      } else {
        triangles.push_back(Triangle(a, b, c));
        triangles.push_back(Triangle(b, c, d));
      }
    }
  }
  return triangles;
}

// Question 44:
// Retourne le produit vectoriel W = A_eta * V.
vector<double> matvec(vector<double> v, vector<Triangle> triangles,
                      Maillage maille) {
  int N = maille.getN(), M = maille.getM();
  VApprox V = VApprox(v, maille);
  V.extendVec();
  vector<double> vv = V.getvGlb();
  vector<double> ww((N + 1) * (M + 1), 0);
  for (Triangle triangle : triangles) {
    vector<vector<double>> BT = triangle.CalcMatBT();
    vector<Noeud> noeuds = triangle.getNoeuds();
    for (int i = 0; i < 3; i++) {
      Noeud ns = noeuds[i];
      int s = ns.numgb(maille);
      double res = 0;
      for (int j = 0; j < 3; j++) {
        Noeud nr = noeuds[j];
        int r = nr.numgb(maille);
        double prod2 = epsilon * triangle.DiffTerm()[j][i] +
                       gama * triangle.ConvectTerm()[j][i] +
                       lambda * triangle.ReacTerm()[j][i];
        res += vv[r] * prod2;
      }
      ww[s] += res;
    }
  }
  V.setvGlb(ww);
  V.IntVec();
  return V.getvInt();
}

// Question 45:
// Retourne le second membre B_eta du système linéaire A_eta * X = B_eta.
vector<double> scdmembre(double (*rhfs)(double, double),
                         vector<Triangle> triangles, Maillage maille) {
  int N = maille.getN(), M = maille.getM();
  vector<double> B((N + 1) * (M + 1), 0);
  for (Triangle triangle : triangles) {
    vector<Noeud> noeuds = triangle.getNoeuds();
    for (int i = 0; i < 3; i++) {
      if (!noeuds[i].est_sur_le_bord(maille)) {
        vector<vector<double>> BT =
            triangle.CalcMatBT({i, (i + 1) % 3, (i + 2) % 3});
        double res = 0;
        vector<double> FT = {BT[0][0] / 2 + noeuds[i].getx(),
                             BT[1][0] / 2 + noeuds[i].gety()};
        res += rhfs(FT[0], FT[1]) / 12;
        FT = {BT[0][1] / 2 + noeuds[i].getx(), BT[1][1] / 2 + noeuds[i].gety()};
        res += rhfs(FT[0], FT[1]) / 12;
        B[noeuds[i].numint(maille)] += res;
      }
    }
  }
  return B;
}

// Fonction pour les tests.
double rhfs_test(double x, double y) {
  return 2 * x * x * x + 3 * y * y * y + 5 * x * y * y + 1;
}

int main(void) {
  // Un petit test pour les noeuds.
  int N = 5;
  int M = 4;
  // int I = (N - 1) * (M - 1);
  // int G = (N + 1) * (M + 1);
  Maillage maille(N, M, 10, 10);
  Noeud noeud(1, 2, maille);
  cout << noeud.numgb(maille) << endl;
  cout << noeud.numint(maille) << endl;
  cout << noeud.num_gb_int(maille, 21) << endl;
  cout << noeud.num_int_gb(maille, 6) << endl;

  cout << endl;

  // Un petit test pour la triangulation.
  vector<Triangle> triangulation = maillageTR(maille);
  cout << "Test de la triangulation pour N = " << N << " et M = " << M << endl;
  for (Triangle triangle : triangulation)
    triangle.affiche_sommets_glb(maille);
  cout << "Il y a " << triangulation.size() << " triangles." << endl;

  cout << endl;
  Triangle triangle(Noeud(0, 0, maille), Noeud(1, 0, maille),
                    Noeud(0, 1, maille));
  cout << triangle.DetMatBT() << endl;

  // Un petit test pour VApprox.
  vector<double> vInt = {1, 3, 2, 4, 5, 9, 0, 8, 6, 7, 4, 5};
  VApprox vec = VApprox(vInt, maille);
  vec.extendVec();
  vector<double> vGlb = vec.getvGlb();
  for (int j = M; j >= 0; j--) {
    for (int i = 0; i <= N; i++) {
      cout << vGlb[(N + 1) * j + i] << " ";
    }
    cout << endl;
  }
  cout << endl;
  vec.IntVec();
  for (int i = 0; i < (N - 1) * (M - 1); i++) {
    if (vInt[i] != vec.getvInt()[i]) {
      cout << "erreur IntVec, v.vInt[" << i << "] = " << vec.getvInt()[i]
           << " , vInt[" << i << "] = " << vInt[i] << endl;
    }
  }

  // Un petit test pour matvec.
  vector<double> vect_test = matvec(vec.getvInt(), triangulation, maille);
  for (double v : vect_test)
    cout << v << " ";
  cout << endl;

  // Un petit test pour scdmembre.
  vector<double> B = scdmembre(rhfs_test, triangulation, maille);
  for (double b : B)
    cout << b << " ";
  cout << endl;

  return 0;
}