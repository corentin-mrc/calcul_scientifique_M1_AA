#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>

#include "vApprox.h"

using namespace std;

VApprox::VApprox(void) {}

VApprox::VApprox(vector<double> vInt, Maillage maille) {
  this->vInt = vInt;
  this->maille = maille;
  this->extendVec();
}

vector<double> VApprox::getvInt(void) { return vInt; }

vector<double> VApprox::getvGlb(void) { return vGlb; }

void VApprox::setvGlb(vector<double> vGlb) { this->vGlb = vGlb; }

void VApprox::extendVec(void) {
  vGlb = vector<double>();
  for (int i = 0; i <= maille.getN(); i++)
    vGlb.push_back(0);
  int numInt = 0;
  for (int j = 1; j < maille.getM(); j++) {
    vGlb.push_back(0);
    for (int i = 1; i < maille.getN(); i++) {
      vGlb.push_back(vInt[numInt]);
      numInt++;
    }
    vGlb.push_back(0);
  }
  for (int i = 0; i <= maille.getN(); i++)
    vGlb.push_back(0);
}

void VApprox::IntVec(void) {
  vInt = vector<double>();
  for (int j = 1; j < maille.getM(); j++) {
    for (int i = 1; i < maille.getN(); i++) {
      vInt.push_back(vGlb[(maille.getN() + 1) * j + i]);
    }
  }
}

vector<double> VApprox::IntVec(vector<double> vGlb, Maillage maille) {
  vector<double> vInt;
  for (int j = 1; j < maille.getM(); j++) {
    for (int i = 1; i < maille.getN(); i++) {
      vInt.push_back(vGlb[(maille.getN() + 1) * j + i]);
    }
  }
  return vInt;
}

double VApprox::normL2(vector<Triangle> triangles) {
  int N = maille.getN(), M = maille.getM();
  this->extendVec();
  vector<double> ww((N + 1) * (M + 1), 0);
  for (Triangle triangle : triangles) {
    vector<vector<double>> BT = triangle.CalcMatBT();
    vector<Noeud> noeuds = triangle.getNoeuds();
    for (int i = 0; i < 3; i++) {
      int s = noeuds[i].numgb(maille);
      double res = 0;
      for (int j = 0; j < 3; j++) {
		int r = noeuds[j].numgb(maille);
        res += vGlb[r] * triangle.ReacTerm()[j][i];
	  }
      ww[s] += res;
    }
  }
  vector<double> AV = VApprox::IntVec(ww, maille);
  // produit scalaire de vInt et AV
  return inner_product(vInt.begin(), vInt.end(), AV.begin(), 0.0);
}

double VApprox::normL2Grad(vector<Triangle> triangles) {
  int N = maille.getN(), M = maille.getM();
  this->extendVec();
  vector<double> ww((N + 1) * (M + 1), 0);
  for (Triangle triangle : triangles) {
    vector<vector<double>> BT = triangle.CalcMatBT();
    vector<Noeud> noeuds = triangle.getNoeuds();
    for (int i = 0; i < 3; i++) {
      int s = noeuds[i].numgb(maille);
      double res = 0;
      for (int j = 0; j < 3; j++) {
		int r = noeuds[j].numgb(maille);
        res += vGlb[r] * triangle.DiffTerm()[j][i];
	  }
      ww[s] += res;
    }
  }
  vector<double> AV = VApprox::IntVec(ww, maille);
  // produit scalaire de vInt et AV
  return inner_product(vInt.begin(), vInt.end(), AV.begin(), 0.0);
}