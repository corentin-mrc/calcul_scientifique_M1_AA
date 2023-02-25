from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, exp, sin

a = 1
b = 1
N = 6
M = 6
divx = 50
divy = 50
epsilon = 1
gama = 1
lambd = 1
I = (N - 1) * (M - 1)
h1 = (2 * a) / N
h2 = (2 * b) / M
p1 = N / (2 * a)
p2 = M / (2 * b)

CA = (gama / epsilon - sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
CB = (gama / epsilon + sqrt(gama**2 / (epsilon**2) + 4 * pi**2 + 4 * lambd / epsilon)) / 2
C1 = 1 / (exp(-CA) - exp(CA - 2 * CB))
C2 = 1 / (exp(-CB) - exp(CB - 2 * CA))

A = np.zeros((I, I))
A_eta_comp = [dict() for i in range(I)]
B = np.zeros(I)
app = [0, 0, 0, 0, 0, 0]

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

def sur_le_bord(i, j):
	return i == 0 or i == N or j == 0 or j == M


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

def ajout(k, j, val):
	if j in A_eta_comp[k]:
		A_eta_comp[k][j] += val
	else:
		A_eta_comp[k][j] = val

def calc_A_comp(i0, j0, i1, j1, i2, j2):
	M_e = [[p1**2 + p2**2, -p1**2, -p2**2], [-p1**2, p1**2, 0], [-p2**2, 0, p2**2]]
	k0 = num_int(i0, j0)
	k1 = num_int(i1, j1)
	k2 = num_int(i2, j2)
	iss = [i0, i1, i2]
	jss = [j0, j1, j2]
	kss = [k0, k1, k2]
	gradx = [(i0 - i1) * p1, (i1 - i0) * p1, 0]
	for k in range(3):
		if not sur_le_bord(iss[k], jss[k]):
			for j in range(3):
				if not sur_le_bord(iss[j], jss[j]):
					aaj = epsilon * M_e[k][j] * h1 * h2 / 2
					aaj += gama * gradx[j] * h1 * h2 / 6
					aaj += lambd * h1 * h2 / (12 if (k == j) else 24)
					if kss[k] == 0 and kss[j] == 0:
						app[3] += epsilon * M_e[k][j] * h1 * h2 / 2
						print("appem :", kss[k], kss[j], epsilon * M_e[k][j] * h1 * h2 / 2)
						app[4] += gama * gradx[j] * h1 * h2 / 6
						app[5] += lambd * h1 * h2 / (12 if (k == j) else 24)
					ajout(kss[k], kss[j], aaj)

def matvec(Xv):
	Yv = np.zeros(I)
	for k in range(I):
		for j in A_eta_comp[k]:
			Yv[k] += A_eta_comp[k][j] * Xv[j]
	return Yv

def inv_syst():
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

def inv_syst_mat_vec():
	X0 = np.ones(I)
	R0 = B - matvec(X0)
	R0_et = R0
	W0 = R0
	for j in range(20):
		AW0 = matvec(W0)
		alpha0 = R0.dot(R0_et) / AW0.dot(R0_et)
		S0 = R0 - alpha0 * AW0
		AS0 = matvec(S0)
		omega0 = AS0.dot(S0) / AS0.dot(AS0)
		X1 = X0 + alpha0 * W0 + omega0 * S0
		R1 = S0 - omega0 * AS0
		beta0 = R1.dot(R0_et) / R0.dot(R0_et) * alpha0 / omega0
		W1 = R1 + beta0 * (W0 - omega0 * AW0)
		R0 = R1
		W0 = W1
		X0 = X1
	return X0


# A_eta, la matrice du système:
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
					if k == 0 and j == 0:
						app[0] += (epsilon * w_j_grad[0] * w_k_grad[0]) * ((h1 / divx) * (h2 / divy))
						app[0] += (epsilon * w_j_grad[1] * w_k_grad[1]) * ((h1 / divx) * (h2 / divy))
						app[1] += (gama * w_j_grad[0] * w_k) * ((h1 / divx) * (h2 / divy))
						app[2] += (lambd * w_j * w_k) * ((h1 / divx) * (h2 / divy))
					y += h2 / divy
				x += h1 / divx
			A[k][j] *= ((h1 / divx) * (h2 / divy))

# A_eta_comp, la matrice du système compressee:
for i in range(N):
	for j in range(M):
		if (i + j) % 2 == 0:
			# C
			i0, j0 = i + 1, j
			i1, j1 = i, j
			i2, j2 = i + 1, j + 1
			calc_A_comp(i0, j0, i1, j1, i2, j2)
			# D
			i0, j0 = i, j + 1
			i1, j1 = i + 1, j + 1
			i2, j2 = i, j
			calc_A_comp(i0, j0, i1, j1, i2, j2)
		else:
			# A
			i0, j0 = i, j
			i1, j1 = i + 1, j
			i2, j2 = i, j + 1
			calc_A_comp(i0, j0, i1, j1, i2, j2)
			# B
			i0, j0 = i + 1, j + 1
			i1, j1 = i, j + 1
			i2, j2 = i + 1, j
			calc_A_comp(i0, j0, i1, j1, i2, j2)

# B_eta, le second membre:
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
			B[k] += (feta(x, y) * wk(x, y, k)) * ((h1 / divx) * (h2 / divy))
			y += h2 / divy
		x += h1 / divx

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
print("\nla matrice A (A_eta) :")
for i in range(19):
	for j in range(19):
		print("{:10.4f}".format(A[i][j]), end = "")
	print()
print("\nA_eta_comp :")
for k in range(19):
	print(k, end = " : " if k > 9 else "  : ")
	for j in A_eta_comp[k]:
		print(j, ("{:10.4f}".format(A_eta_comp[k][j])), end = " | ")
	print()
print("\nA_eta * w_eta")
print((A.dot(w_exa))[:30])
print("\nA_eta * X0 :")
print(AX0[:30])
print("\nA_eta * B_eta :")
print(AB_eta[:30])
w_appr = inv_syst_mat_vec()
print("\nw solution approchee mat vec :")
print(w_appr[:30])
w_appr = inv_syst()
print("\nw solution approchee :")
print(w_appr[:30])
w_numpy = np.linalg.solve(A, B)
print("\nw solution par numpy.linalg.solve :")
print(w_numpy[:30])
print("\nA[0][0] = ")
print(app[0], app[1], app[2], app[0] + app[1] + app[2])
print(app[3], app[4], app[5], app[3] + app[4] + app[5])


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

