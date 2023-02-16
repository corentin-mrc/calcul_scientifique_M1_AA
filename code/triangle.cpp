#include <cmath>
#include <iostream>
#include <vector>

#include "triangle.h"

using namespace std;

Triangle::Triangle(void) {}

Triangle::Triangle(Noeud n0, Noeud n1, Noeud n2) {
  this->noeuds = {n0, n1, n2};
}

Triangle::Triangle(vector<Noeud> noeuds) { this->noeuds = noeuds; }

vector<Noeud> Triangle::get_noeuds(void) { return noeuds; }

vector<vector<double>> Triangle::calc_mat_BT() {
  return {{noeuds[1].get_x() - noeuds[0].get_x(),
           noeuds[2].get_x() - noeuds[0].get_x()},
          {noeuds[1].get_y() - noeuds[0].get_y(),
           noeuds[2].get_y() - noeuds[0].get_y()}};
}

vector<vector<double>> Triangle::calc_mat_BT(vector<int> permut) {
  return {{noeuds[permut[1]].get_x() - noeuds[permut[0]].get_x(),
           noeuds[permut[2]].get_x() - noeuds[permut[0]].get_x()},
          {noeuds[permut[1]].get_y() - noeuds[permut[0]].get_y(),
           noeuds[permut[2]].get_y() - noeuds[permut[0]].get_y()}};
}

double Triangle::det_mat_BT(void) {
  vector<vector<double>> BT = calc_mat_BT();
  return BT[0][0] * BT[1][1] - BT[0][1] * BT[1][0];
}

vector<vector<double>> Triangle::inv_mat_BT(void) {
  vector<vector<double>> BT = calc_mat_BT();
  double det = det_mat_BT();
  return {{BT[1][1] / det, -BT[0][1] / det}, {-BT[1][0] / det, BT[0][0] / det}};
}

vector<vector<double>> Triangle::diff_terme(void) {
  vector<vector<double>> BI = inv_mat_BT();
  double d = abs(det_mat_BT()) / 2;
  double d1 = BI[0][0] + BI[1][0];
  double d2 = BI[0][1] + BI[1][1];
  double m00 = (d1 * d1 + d2 * d2) * d;
  double m01 = (-BI[0][0] * d1 - BI[0][1] * d2) * d;
  double m02 = (-BI[1][0] * d1 - BI[1][1] * d2) * d;
  double m11 = (BI[0][0] * BI[0][0] + BI[0][1] * BI[0][1]) * d;
  double m12 = (BI[0][0] * BI[1][0] + BI[0][1] * BI[1][1]) * d;
  double m22 = (BI[1][0] * BI[1][0] + BI[1][1] * BI[1][1]) * d;
  return {{m00, m01, m02}, {m01, m11, m12}, {m02, m12, m22}};
}

vector<vector<double>> Triangle::convect_terme(void) {
  double c = abs(det_mat_BT()) / 6;
  return {{-c, -c, -c}, {c, c, c}, {0, 0, 0}};
}

vector<vector<double>> Triangle::react_terme(void) {
  double det = abs(det_mat_BT());
  double r1 = det / 12;
  double r2 = det / 24;
  return {{r1, r2, r2}, {r2, r1, r2}, {r2, r2, r1}};
}