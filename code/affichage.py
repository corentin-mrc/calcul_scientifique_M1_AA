from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import floor, sqrt, pi, exp, sin


N = 25
epsilon = 1
M = N
I = (N - 1)**2
a = 1
b = a
h1 = 2 * a / N
h2 = 2 * b / M
gama = 1
lambd = 1
det = h1 * h2

A = (gama / epsilon - sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
B = (gama / epsilon + sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
C1 = 1 / (exp(-A) - exp(A - 2 * B))
C2 = 1 / (exp(-B) - exp(B - 2 * A))

def solution_exacte(x, y):
	return (C1 * np.exp(A * x) + C2 * np.exp(B * x)) * np.sin(pi * y)

def u_star(x, y):
	return (a - x) * np.sin(pi * y) / (2 * a)


filename = "solution_approchee.txt"
with open(filename) as f:
    content = f.readlines()
    wk = []
    for line in content:
        if line != '' and line != '\n':
            wk.append(float(line))
# w_eta_h de i et j
fxy = [[0] * (M + 1) for i in range(N + 1)]
for k in range(I):
	i = (k % (N - 1)) + 1
	j = (k // (N - 1)) + 1
	xi = a * (2 * i - N) / N
	yj = b * (2 * j - M) / M
	fxy[i][j] = wk[k] + u_star(xi, yj)
X_a = [[0] * (M + 1) for i in range(N + 1)]
Y_a = [[0] * (M + 1) for i in range(N + 1)]
for i in range(N + 1):
	xi = a * (2 * i - N) / N
	for j in range(M + 1):
		yj = b * (2 * j - M) / M
		X_a[i][j] = xi
		Y_a[i][j] = yj

def f(x, y):
	i = floor(N * (x + a) / (2 * a))
	xi = (2 * i - N) * a / N
	j = floor(M * (y + b) / (2 * b))
	yj = (2 * j - M) * b / M
	if (i + j) % 2 == 0:
		i0 = i; j0 = j; x0 = xi; y0 = yj
		if (x - x0) / h1 < (y - y0) / h2:
			i1 = i0; j1 = j0 + 1; x1 = x0; y1 = y0 + h2
		else:
			i1 = i0 + 1; j1 = j0; x1 = x0 + h1; y1 = y0
		i2 = i0 + 1; j2 = j0 + 1; x2 = x0 + h1; y2 = y0 + h2
	else:
		i0 = i + 1; j0 = j; x0 = xi + h1; y0 = yj
		if (x0 - x) / h1 < (y - y0) / h2:
			i1 = i0; j1 = j0 + 1; x1 = x0; y1 = y0 + h2
		else:
			i1 = i0 - 1; j1 = j0; x1 = x0 - h1; y1 = y0
		i2 = i0 - 1; j2 = j0 + 1; x2 = x0 - h1; y2 = y0 + h2
	#Bt = [[x1 - x0, x2 - x0], [y1 - y0, y2 - y0]]
	Bi = [[(y2 - y0) / det, (x0 - x2) / det], [(y0 - y1) / det, (x1 - x0) / det]]
	xchap = Bi[0][0] * (x - x0) + Bi[0][1] * (y - y0)
	ychap = Bi[1][0] * (x - x0) - Bi[1][1] * (y - y0)
	lambdas = [1 - xchap - ychap, xchap, ychap]
	return lambdas[0] * fxy[i0][j0] + lambdas[1] * fxy[i1][j1] + lambdas[2] * fxy[i2][j2]


fig = plt.figure()

nb_pt = 60
xp = np.linspace(-a, a, nb_pt)
yp = np.linspace(-b, b, nb_pt)

X, Y = np.meshgrid(xp, yp)
Z = np.zeros((nb_pt, nb_pt))
for i in range(nb_pt):
	x = a * (2 * i - nb_pt) / nb_pt
	for j in range(nb_pt):
		y = b * (2 * j - nb_pt) / nb_pt
#		Z[i][j] = f(x, y) + u_star(x, y)
#		Z[i][j] = u_star(x, y)
ax = plt.axes(projection='3d')
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_surface(np.array(X_a), np.array(Y_a), np.array(fxy), rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_title('solution exacte')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x, y)')


plt.show()


