#include <iostream>
#include <vector>
#include "constantes.h"

using namespace std;

// question 3
double fsecond_membre(double (*ug) (double), double (*ud) (double), double (*ugpp) (double), double (*udpp) (double), double x, double y)
{
	return (-epsilon * (a - x) * ugpp(y) - epsilon * (a + x) * udpp(y) - gama * ug(y) + gama * ud(y) + lambda * (a - x) * ug(y) + lambda * (a + x) * ud(y)) / (2 * a);
}

// question 9
vector<double> subdiv(int N)
{
	vector<double> subdiv = vector<double>(N + 1);
	for (int i = 0; i <= N; i++)
	{
		subdiv[i] = -a + (2.0 * i * a) / N;
	}
	return subdiv;
}

// question 11
int numgb(int N, int M, int i, int j)
{
	return (M + 1) * i + j;
}

// question 13
void invnumgb(int N, int M, int s, int& i, int& j)
{
	i = s / (M + 1);
	j = s % (M + 1);
}

int main()
{
	cout << "Hello World! " << a << " , " << b << " , " << epsilon << " , " << gama << " , " << lambda << endl;
	int N = 15;
	vector<double> test9 = subdiv(N);
	for (int i = 0; i <= N; i++)
		cout << test9[i] << " ";
	cout << endl;
	return 0;
}