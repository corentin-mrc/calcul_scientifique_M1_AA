from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, exp, sin

a = 1
b = 1
N = 12
M = 12
divx = 50
divy = 50
epsilon = 1
gama = 1
lambd = 1
h1 = (2 * a) / N
h2 = (2 * b) / M
I = (N - 1) * (M - 1)

CA = (gama / epsilon - sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
CB = (gama / epsilon + sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
C1 = 1 / (exp(-CA) - exp(CA - 2 * CB))
C2 = 1 / (exp(-CB) - exp(CB - 2 * CA))

def u_exact(x, y):
	return (C1 * exp(CA * x) + C2 * exp(CB * x)) * sin(pi * y)

def u_star(x, y):
	return (a - x) * sin(pi * y) / (2 * a)

def w_exact(x, y):
	return u_exact(x, y) - u_star(x, y)


def feta(x, y):
	return ((pi**2 * epsilon + lambd) * x + gama - lambd * a - epsilon * a * pi**2) * sin(pi * y) / (2 * a)

def num_int(i, j):
	return (j - 1) * (N - 1) + (i - 1)

def inv_num_int(k):
	return (k % (N - 1)) + 1, (k // (N - 1)) + 1


def int_coord(k):
	i, j = inv_num_int(k)
	return (2 * i - N) * a / N, (2 * j - M) * b / M


def wk(x, y, k):
	i, j = (k % (N - 1)) + 1, (k // (N - 1)) + 1
	xi, yj = (2 * i - N) * a / N, (2 * j - M) * b / M
	if (i + j) % 2 == 0:
		if abs(xi - x) < h1 and abs(yj - y) < h2:
			return 1 - max(abs(xi - x) / h1, abs(yj - y) / h2)
		return 0
	if abs(xi - x) / h1 + abs(yj - y) / h2 < 1:
		return 1 - abs(xi - x) / h1 - abs(yj - y) / h2
	return 0

def grad_wk(x, y, k):
	i, j = (k % (N - 1)) + 1, (k // (N - 1)) + 1
	xi, yj = (2 * i - N) * a / N, (2 * j - M) * b / M
	if (i + j) % 2 == 0:
		if abs(xi - x) < h1 and abs(yj - y) < h2:
			xa = abs(x - xi) / h1
			yb = abs(y - yj) / h2
			return 0 if (yb > xa) else (-1 / h1 if x > xi else 1 / h1), 0 if (xa > yb) else (-1 / h2 if y > yj else 1 / h2)
		return 0, 0
	if abs(xi - x) / h1 + abs(yj - y) / h2 < 1:
		return -1 / h1 if (x > xi) else 1 / h1, -1 / h2 if (y > yj) else 1 / h2
	return 0, 0

# A_eta, la matrice du système:
A = np.zeros((I, I))
for k in range(I):
	ik, jk = inv_num_int(k)
	xik, yjk = int_coord(k)
	for ij in range(max(1, ik - 1), min(N - 1, ik + 1) + 1):
		for jj in range(max(1, jk - 1), min(M - 1, jk + 1) + 1):
			j = num_int(ij, jj)
			x = xik - h1
			while x < xik + h1:
				y = yjk - h2
				while y < yjk + h2:
					w_k_grad = grad_wk(x, y, k)
					w_k = wk(x, y, k)
					w_j_grad = grad_wk(x, y, j)
					w_j = wk(x, y, j)
					A[k][j] += epsilon * w_j_grad[0] * w_k_grad[0]
					A[k][j] += epsilon * w_j_grad[1] * w_k_grad[1]
					A[k][j] += gama * w_j_grad[0] * w_k
					A[k][j] += lambd * w_j * w_k
					y += h2 / divy
				x += h1 / divx
			A[k][j] /= divx * divy

# B_eta, le second membre:
B = np.zeros(I)
for k in range(I):
	apr_inte = 0
	i = (k % (N - 1)) + 1
	j = (k // (N - 1)) + 1
	xi = (2 * i - N) * a / N
	yj = (2 * j - M) * b / M
	x = xi - h1
	while x < xi + h1:
		y = yj - h2
		while y < yj + h2:
			B[k] += feta(x, y) * wk(x, y, k) / divx / divy
			y += h2 / divy
		x += h1 / divx

def inv_syst(A, B):
	X0 = np.ones(I)
	R0 = B - A.dot(X0)
	R0_et = R0
	W0 = R0
	for j in range(20):
		AW0 = A.dot(W0)
		alpha0 = R0.dot(R0_et) / AW0.dot(R0_et)
		S0 = R0 - alpha0 * AW0
		AS0 = A.dot(S0)
		omega0 = AS0.dot(S0) / AS0.dot(AS0)
		X1 = X0 + alpha0 * W0 + omega0 * S0
		R1 = S0 - omega0 * AS0
		beta0 = R1.dot(R0_et) / R0.dot(R0_et) * alpha0 / omega0
		W1 = R1 + beta0 * (W0 - omega0 * AW0)
		R0 = R1
		W0 = W1
		X0 = X1
	return X0

X0 = [1] * I
AX0 = A.dot(X0)
AB_eta = A.dot(B)
w_exa = np.zeros(I)
for k in range(I):
	xi, yj = int_coord(k)
	w_exa[k] = w_exact(xi, yj)
print("A_eta :")
print(A[:10][:10])
print("\nle w exact :")
print(w_exa[:30])
print("\nsecond membre (B_eta) :")
print(B[:30])
print("\nA_eta * w_eta")
for i in range(19):
	for j in range(19):
		print("{:10.4f}".format(A[i][j]), end = "")
	print()
print((A.dot(w_exa))[:30])
print("\nA_eta * X0 :")
print(AX0[:30])
print("\nA_eta * B_eta :")
print(AB_eta[:30])
w_appr = inv_syst(A, B)
print("\nw solution approchee :")
print(w_appr[:30])
w_numpy = np.linalg.solve(A, B)
print("\nw solution par numpy.linalg.solve :")
print(w_numpy[:30])

X = np.zeros((N + 1, M + 1))
Y = np.zeros((N + 1, M + 1))
Z_w_appr = np.zeros((N + 1, M + 1))
Z_u_appr = np.zeros((N + 1, M + 1))
Z_w_exa = np.zeros((N + 1, M + 1))
Z_u_exa = np.zeros((N + 1, M + 1))

for i in range(N + 1):
	xi = (2 * i - N) * a / N
	for j in range(M + 1):
		yj = (2 * j - M) * b / M
		X[i][j] = xi
		Y[i][j] = yj
		if i != 0 and i != N and j != 0 and j != M:
			Z_w_appr[i][j] = w_appr[num_int(i, j)]
		Z_u_appr[i][j] = Z_w_appr[i][j] + u_star(xi, yj)
		Z_w_exa[i][j] = w_exact(xi, yj)
		Z_u_exa[i][j] = u_exact(xi, yj)
# choisir quoi afficher entre     0 : w_appr   |   1 : u_appr    |    2 : w_exa    |    3 : u_exa
choix = 0
Z = [Z_w_appr, Z_u_appr, Z_w_exa, Z_u_exa]
Z_nom = ["w_appr", "u_appr", "w_exa", "u_exa"]
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z[choix], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_title(Z_nom[choix])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel(Z_nom[choix] + "(x, y)")

plt.show()

