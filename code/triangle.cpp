#include <cmath>
#include <iostream>
#include <vector>

#include "triangle.h"

using namespace std;

Triangle::Triangle(void) {}

Triangle::Triangle(Noeud n0, Noeud n1, Noeud n2) {
  this->n0 = n0;
  this->n1 = n1;
  this->n2 = n2;
}

vector<vector<double>> Triangle::CalcMatBT(void) {
  return {{n1.getx() - n0.getx(), n2.getx() - n0.getx()},
          {n1.gety() - n0.gety(), n2.gety() - n0.gety()}};
}

double Triangle::DetMatBT(void) {
  vector<vector<double>> BT = CalcMatBT();
  return BT[0][0] * BT[1][1] - BT[0][1] * BT[1][0];
}

vector<vector<double>> Triangle::InvMatBT(void) {
  vector<vector<double>> BT = CalcMatBT();
  double det = DetMatBT();
  return {{BT[1][1] / det, -BT[0][1] / det}, {-BT[1][0] / det, BT[0][0] / det}};
}

vector<vector<double>> Triangle::DiffTerm(void) {
  vector<vector<double>> BI = InvMatBT();
  double d = DetMatBT() / 2;
  //double d1 = pow(BI[0][0] + BI[0][1], 2) * d;
  //double d2 = pow(BI[1][0] + BI[1][1], 2) * d;
  //return {{d1 + d2, -d1, -d2}, {-d1, d1, 0}, {-d2, 0, d2}};
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

vector<vector<double>> Triangle::ConvectTerm(void) {
  double c = abs(DetMatBT()) / 6;
  return {{-c, -c, -c}, {c, c, c}, {0, 0, 0}};
}

vector<vector<double>> Triangle::ReacTerm(void) {
  double det = abs(DetMatBT());
  double r1 = det / 12;
  double r2 = det / 24;
  return {{r1, r2, r2}, {r2, r1, r1}, {r2, r1, r1}};
}

void Triangle::affiche_sommets_glb(Maillage maille) {
  cout << n0.numgb(maille) << " " << n1.numgb(maille) << " " << n2.numgb(maille)
       << endl;
}