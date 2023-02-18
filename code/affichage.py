from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import floor, sqrt, pi, exp, sin


epsilon = 1
a = 1
b = 1
gama = 1
lambd = 1

A = (gama / epsilon - sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
B = (gama / epsilon + sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
C1 = 1 / (exp(-A) - exp(A - 2 * B))
C2 = 1 / (exp(-B) - exp(B - 2 * A))

def solution_exacte(x, y):
	return (C1 * np.exp(A * x) + C2 * np.exp(B * x)) * np.sin(pi * y)

def u_star(x, y):
	return (a - x) * np.sin(pi * y) / (2 * a)


#filename = "solution_exacte.txt"
filename = "solution_approchee.txt"
with open(filename) as f:
    content = f.readlines()
    wk = []
    for line in content:
        if line != '' and line != '\n':
            wk.append(float(line))
I = len(wk)
# on sppose M et N egaux
N = round(sqrt(I) + 1)
M = N
h1 = 2 * a / N
h2 = 2 * b / M
det = h1 * h2
# w_eta_h de i et j
fxy = [[0] * (M + 1) for i in range(N + 1)]
for k in range(I):
	i = (k % (N - 1)) + 1
	j = (k // (N - 1)) + 1
	xi = a * (2 * i - N) / N
	yj = b * (2 * j - M) / M
	fxy[i][j] = wk[k]# + u_star(xi, yj)
#for j in range(M + 1):
#	yj = b * (2 * j - M) / M
#	fxy[0][j] = u_star(-a, yj)
X_a = [[0] * (M + 1) for i in range(N + 1)]
Y_a = [[0] * (M + 1) for i in range(N + 1)]
for i in range(N + 1):
	xi = a * (2 * i - N) / N
	for j in range(M + 1):
		yj = b * (2 * j - M) / M
		X_a[i][j] = xi
		Y_a[i][j] = yj
#print(fxy)
#fxy = [[0, 0, 0], [0, .5, 0], [0, 0, 0]]

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
	det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
	Bi = [[(y2 - y0) / det, (x0 - x2) / det], [(y0 - y1) / det, (x1 - x0) / det]]
	xchap = Bi[0][0] * (x - x0) + Bi[0][1] * (y - y0)
	ychap = Bi[1][0] * (x - x0) + Bi[1][1] * (y - y0)
	lambdas = [1 - xchap - ychap, xchap, ychap]
	return lambdas[0] * fxy[i0][j0] + lambdas[1] * fxy[i1][j1] + lambdas[2] * fxy[i2][j2]

def f2(x, y):
	i = round(N * (x + a) / (2 * a))
	xi = (2 * i - N) * a / N
	if xi > x:
		i -= 1
		xi -= h1
	j = round(M * (y + b) / (2 * b))
	yj = (2 * j - M) * b / M
	if yj > y:
		j -= 1
		yj -= h2
	i0 = 0
	i1 = 0
	i2 = 0
	j0 = 0
	j1 = 0
	j2 = 0
	x0 = .1
	x1 = .1
	x2 = .1
	y0 = .1
	y1 = .1
	y2 = .1
	if (((i + j) % 2) == 0):
		i0 = i
		j0 = j
		x0 = xi
		y0 = yj
		if (((x - x0) / h1) < ((y - y0) / h2)):
			i1 = i0
			j1 = j0 + 1
			x1 = x0
			y1 = y0 + h2
		else:
			i1 = i0 + 1
			j1 = j0
			x1 = x0 + h1
			y1 = y0
		i2 = i0 + 1
		j2 = j0 + 1
		x2 = x0 + h1
		y2 = y0 + h2
	else:
		i0 = i + 1
		j0 = j
		x0 = xi + h1
		y0 = yj
		if ((x0 - x) / h1) < ((y - y0) / h2):
			i1 = i0
			j1 = j0 + 1
			x1 = x0
			y1 = y0 + h2
		else:
			i1 = i0 - 1
			j1 = j0
			x1 = x0 - h1
			y1 = y0
		i2 = i0 - 1
		j2 = j0 + 1
		x2 = x0 - h1
		y2 = y0 + h2
	#Bt = [[x1 - x0, x2 - x0], [y1 - y0, y2 - y0]]
	Bi = [[(y2 - y0) / det, (x0 - x2) / det], [(y0 - y1) / det, (x1 - x0) / det]]
	xchap = (Bi[0][0] * (x - x0)) + (Bi[0][1] * (y - y0))
	ychap = (Bi[1][0] * (x - x0)) + (Bi[1][1] * (y - y0))
	lambdas = [1 - xchap - ychap, xchap, ychap]
	lambda0 = (y1 - y2) * x + (x2 - x1) * y + x1 * y2 - x2 * y1
	lambda1 = (y2 - y0) * x + (x0 - x2) * y + x2 * y0 - x0 * y2
	lambda2 = (y0 - y1) * x + (x1 - x0) * y + x0 * y1 - x1 * y0
	dete = x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2
	lambdas2 = [lambda0 / dete, lambda1 / dete, lambda2 / dete]
	err = ""
	if abs(lambdas[0] - lambdas2[0]) > .001 or abs(lambdas[1] - lambdas2[1]) > .001 or abs(lambdas[1] - lambdas2[1]) > .001:
		err += "  lambdas differents"
	if lambdas[0] > 1.001 or lambdas[1] > 1.001 or lambdas[2] > 1.001:
		err += "  lambda > 1"
	if err != "" and err == "c":
		print("determinant :", dete)
		print("lambdas1 :", lambdas)
		print("lambdas2 :", lambdas2)
		print("x , y  :", x, "\t", y)
		print("i , j  :", i, "\t", j)
		print("xi, yj :", xi, "\t", yj)
		print("x0, y0 :", x0, "\t", y0)
		print("x1, y1 :", x1, "\t", y1)
		print("x2, y2 :", x2, "\t", y2)
		print("i0, j0 :", i0, "\t", j0)
		print("i1, j1 :", i1, "\t", j1)
		print("i2, j2 :", i2, "\t", j2)
		print(err)
		print()
	return ((lambdas2[0] * fxy[i0][j0]) + (lambdas2[1] * fxy[i1][j1]) + (lambdas2[2] * fxy[i2][j2]))


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
		Z[i][j] = f(x, y)# + u_star(x, y)
#		Z[i][j] = u_star(x, y)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#ax.plot_surface(np.array(X_a), np.array(Y_a), np.array(fxy), rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_title('solution exacte')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x, y)')


plt.show()


