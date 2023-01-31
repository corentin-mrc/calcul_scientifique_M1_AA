#include <cmath>
#include <iostream>
#include <vector>

#include "vApprox.h"

using namespace std;

VApprox::VApprox(void) {}

VApprox::VApprox(vector<double> vInt, Maillage maille) {
  this->vInt = vInt;
  this->maille = maille;
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
