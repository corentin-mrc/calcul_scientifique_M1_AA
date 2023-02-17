#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "boite_a_outils.h"

using namespace std;

vector<double> operator+(vector<double> A, vector<double> B) {
  vector<double> C;
  for (size_t i = 0; i < B.size(); ++i)
    C.push_back(A[i] + B[i]);
  return C;
}

vector<double> operator-(vector<double> A, vector<double> B) {
  vector<double> C;
  for (size_t i = 0; i < B.size(); ++i)
    C.push_back(A[i] - B[i]);
  return C;
}

vector<double> operator*(double scalaire, vector<double> B) {
  vector<double> C;
  for (size_t i = 0; i < B.size(); ++i)
    C.push_back(scalaire * B[i]);
  return C;
}

double operator*(vector<double> A, vector<double> B) {
  double C;
  for (size_t i = 0; i < B.size(); ++i)
    C += A[i] * B[i];
  return C;
}

double max(vector<double> A) {
  double max = 0;
  for (double i : A)
    max = i > 0 ? i > max ? i : max : -i > max ? -i : max;
  return max;
}

double f_second_membre(double (*u_g)(double), double (*u_d)(double),
                       double (*u_gpp)(double), double (*u_dpp)(double),
                       double x, double y) {
  return (epsilon * (a - x) * u_gpp(y) + epsilon * (a + x) * u_dpp(y) +
          gama * u_g(y) - gama * u_d(y) - lambda * (a - x) * u_g(y) -
          lambda * (a + x) * u_d(y)) /
         (2 * a);
}

vector<double> extend_vec(vector<double> V_int) {
  vector<double> V_glb;
  for (int i = 0; i <= N; i++)
    V_glb.push_back(0);
  int numInt = 0;
  for (int j = 1; j < M; j++) {
    V_glb.push_back(0);
    for (int i = 1; i < N; i++) {
      V_glb.push_back(V_int[numInt]);
      numInt++;
    }
    V_glb.push_back(0);
  }
  for (int i = 0; i <= N; i++)
    V_glb.push_back(0);
  return V_glb;
}

vector<double> int_vec(vector<double> V_glb) {
  vector<double> V_int;
  for (int j = 1; j < M; j++) {
    for (int i = 1; i < N; i++)
      V_int.push_back(V_glb[(N + 1) * j + i]);
  }
  return V_int;
}

double norme_L2(vector<double> V, Maillage maille) {
  vector<double> V_glb = extend_vec(V);
  vector<double> WW((N + 1) * (M + 1), 0);
  for (Triangle triangle : maille.get_triangulation()) {
    vector<Noeud> noeuds = triangle.get_noeuds();
    for (int i = 0; i < 3; i++) {
      int s = maille.num_gb_noeud(noeuds[i]);
      double res = 0;
      for (int j = 0; j < 3; j++) {
        int r = maille.num_gb_noeud(noeuds[j]);
        res += V_glb[r] * triangle.react_terme()[j][i];
      }
      WW[s] += res;
    }
  }
  vector<double> AV = int_vec(WW);
  return AV * V;
}

double norme_L2_grad(vector<double> V, Maillage maille) {
  vector<double> V_glb = extend_vec(V);
  vector<double> WW((N + 1) * (M + 1), 0);
  for (Triangle triangle : maille.get_triangulation()) {
    vector<Noeud> noeuds = triangle.get_noeuds();
    for (int i = 0; i < 3; i++) {
      int s = maille.num_gb_noeud(noeuds[i]);
      double res = 0;
      for (int j = 0; j < 3; j++) {
        int r = maille.num_gb_noeud(noeuds[j]);
        res += V_glb[r] * triangle.diff_terme()[j][i];
      }
      WW[s] += res;
    }
  }
  vector<double> AV = int_vec(WW);
  return AV * V;
}

vector<double> mat_vec(vector<double> V, Maillage maille) {
  vector<double> VV = extend_vec(V);
  vector<double> WW((N + 1) * (M + 1), 0);
  for (Triangle triangle : maille.get_triangulation()) {
    vector<Noeud> noeuds = triangle.get_noeuds();
    for (int i = 0; i < 3; i++) {
      int s = maille.num_gb_noeud(noeuds[i]);
      double res = 0;
      for (int j = 0; j < 3; j++) {
        int r = maille.num_gb_noeud(noeuds[j]);
        double prod2 = epsilon * triangle.diff_terme()[j][i] +
                       gama * triangle.convect_terme()[j][i] +
                       lambda * triangle.react_terme()[j][i];
        res += VV[r] * prod2;
      }
      WW[s] += res;
    }
  }
  return int_vec(WW);
}

vector<double> scd_membre(double (*rhfs)(double, double), Maillage maille) {
  vector<double> B((N - 1) * (M - 1), 0);
  for (Triangle triangle : maille.get_triangulation()) {
    vector<Noeud> noeuds = triangle.get_noeuds();
    for (int i = 0; i < 3; i++) {
      if (!maille.est_sur_le_bord(noeuds[i])) {
        vector<vector<double>> BT =
            triangle.calc_mat_BT({i, (i + 1) % 3, (i + 2) % 3});
        double res = 0;
        // FT(1/2, 0) = BT * (1/2, 0) + (x0, y0):
        vector<double> FT = {BT[0][0] / 2 + noeuds[i].get_x(),
                             BT[1][0] / 2 + noeuds[i].get_y()};
        // wk(FT(1/2, 0)) = 1/2
        res += rhfs(FT[0], FT[1]) / 12;
        // FT(0, 1/2) = BT * (0, 1/2) + (x0, y0):
        FT = {BT[0][1] / 2 + noeuds[i].get_x(),
              BT[1][1] / 2 + noeuds[i].get_y()};
        // wk(FT(0, 1/2)) = 1/2
        res += rhfs(FT[0], FT[1]) / 12;
        B[maille.num_int_noeud(noeuds[i])] += res;
      }
    }
  }
  return B;
}

vector<double> inv_syst(vector<double> B_eta, Maillage maille,
                        int max_iteration) {
  vector<double> X0(B_eta.size(), 1);
  vector<double> R0 = B_eta - mat_vec(X0, maille);
  vector<double> R0_etoile = R0;
  vector<double> W0 = R0;
  for (int j = 0; j < max_iteration; j++) {
    vector<double> AW0 = mat_vec(W0, maille);
    double alpha0 = (R0 * R0_etoile) / (AW0 * R0_etoile);
    vector<double> S0 = R0 - (alpha0 * AW0);
    vector<double> AS0 = mat_vec(S0, maille);
    double omega0 = (AS0 * S0) / (AS0 * AS0);
    vector<double> X1 = X0 + (alpha0 * W0) + (omega0 * S0);
    vector<double> R1 = S0 - (omega0 * AS0);
    double beta0 = ((R1 * R0_etoile) / (R0 * R0_etoile)) * (alpha0 / omega0);
    vector<double> W1 = R1 + (beta0 * (W0 - (omega0 * AW0)));
    R0 = R1;
    W0 = W1;
    X0 = X1;
  }
  return X0;
}

vector<double> erreurs(double (*sol_exa)(double, double),
                       vector<double> sol_appr, Maillage maille) {
  vector<double> w;
  int I = (N - 1) * (M - 1);
  for (int k = 0; k < I; k++) {
    vector<double> xy = maille.int_coord(k);
    w.push_back(sol_exa(xy[0], xy[1]));
  }
  vector<double> erreur = w - sol_appr;
  return {norme_L2(erreur, maille) / norme_L2(w, maille),
          norme_L2_grad(erreur, maille) / norme_L2_grad(w, maille),
          max(erreur) / max(w)};
}

double u_eta(double x, double y) {
  double A = (gama / epsilon - sqrt(gama * gama / epsilon / epsilon +
                                    4 * M_PI * M_PI + 4 * lambda / epsilon)) /
             2;
  double B = (gama / epsilon + sqrt(gama * gama / epsilon / epsilon +
                                    4 * M_PI * M_PI + 4 * lambda / epsilon)) /
             2;
  double C_1 = 1 / (exp(-A) - exp(A - 2 * B));
  double C_2 = 1 / (exp(-B) - exp(B - 2 * A));
  double U_x = C_1 * exp(A * x) + C_2 * exp(B * x);
  return U_x * sin(M_PI * y);
}
