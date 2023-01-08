double a = 10;
double b = 5;
double epsilon = 1;
double gama = 1;
double lambda = 1;

class Noeud
{
	int i;
	int j;
	double x;
	double y;
	Noeud(int pi, int pj, double px, double py)
	{
		this->i = pi;
		this->j = pj;
		this->x = px;
		this->y = py;
	}
};