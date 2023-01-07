#include <iostream>
#include <vector>
#include "constantes.h"

using namespace std;

double fsecond_membre(double (*ug) (double), double (*ud) (double), double (*ugpp) (double), double (*udpp) (double), double x, double y)
{
	return (-epsilon * (a - x) * ugpp(y) - epsilon * (a + x) * udpp(y) - gama * ug(y) + gama * ud(y) + lambda * (a - x) * ug(y) + lambda * (a + x) * ud(y)) / (2 * a);
}

int main()
{
	cout << "Hello World! " << a << " , " << b << " , " << epsilon << " , " << gama << " , " << lambda << endl;
	return 0;
}