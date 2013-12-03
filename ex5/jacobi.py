import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import dok_matrix

#construct the A matrix
def build_A(N):
	_A = dok_matrix((N*N,N*N),dtype=sp.float32)
	
	index_array = np.arange(N*N)
	index_array = np.reshape(index_array,(N,N))

	for i in range(N):
		for j in range(N):
			center = index_array[i,j]
			_A[center,center] = -4
			for delta in ((-1,0),(1,0),(0,1),(0,-1)):
				_i = i + delta[0]
				_j = j + delta[1]
				#periodic wrapping
				_i = _i % N
				_j = _j % N
				index = index_array[_i,_j]
				_A[center,index] = 1
	return _A

# This function improves the solution Ax=b with one Jacobi step, with the 
# result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
def jacobi_step(x, b, N):
	_x = x.flatten()
	_b = b.flatten()
	xnew = D.dot(_b) + DLU.dot(_x)
	x = np.reshape(xnew,(N,N))
	return x

# This function calculates the resdiuum vector res = b - Ax, for input vectors
# of length N. The output is stored in res.
def calc_residuum(x, b, N, res):
	_b = b.flatten()
	_x = x.flatten()
	res = _b - A.dot(_x)
	res = np.reshape(res,(N,N))
	return res

# This function calculates the norm of the vector of length N, 
# defined as the usual quadratic vector norm.   
def norm_of_residual(res, N):
    sum = np.sum(res**2)
    return np.sqrt(sum)
    
def iteration_scheme():
	steps = 2000
	L = 1.0
	h = L / N
	eta = 0.1 * L
	rho0 = 10.0

	res = np.zeros((N,N))

	#now set-up the density field
	x = (np.arange(-N/2,N/2)+0.5)*h
	mesh = np.meshgrid(x,x)
	r2 = mesh[0] * mesh[0] + mesh[1] * mesh[1]

	rho = rho0 * np.exp(-r2/(2*eta**2))
	sum = np.sum(rho)

	#initialize the starting values for phi[] and b[]
	rho -= sum/(N**2)
	b = 4*np.pi*h**2*rho
	phi = np.zeros((N,N))

	#open a file for outputting the residuum values, and then do 2000 Jacobi steps
	residuals = np.zeros(steps)
	f = file("res_jacobi.txt", "w")
	f.write("%d\n"%steps)

	for i in np.arange(steps):
		phi = jacobi_step(phi, b, N)

		res = calc_residuum(phi, b, N, res)

		r = norm_of_residual(res, N)

		print("iter=%d:  residual=%g\n"%(i, r))
		f.write("%g\n"%r)
		residuals[i] = r

	f.close()

	plt.semilogy(residuals, label="jacobi")
	plt.show()

#build the A matrix
N = 256
A = build_A(N)
d = A.diagonal()
D = dok_matrix((N*N,N*N),dtype=sp.float32)
D.setdiag(d)
LU = A.copy()
LU = -1.*(LU - D)
d = 1./d
D = dok_matrix((N*N,N*N),dtype=sp.float32)
D.setdiag(d)
DLU = D*LU
iteration_scheme()
