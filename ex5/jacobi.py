import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import dok_matrix

def to_quadrants(a):
		#ul_|_ur  upper left(ul), ect
		#ll | lr
		
		b = N*N/2
		ur = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
		lr = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
		ul = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
		ll = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
		
		#arg!!!
		for entry in a.iteritems():
			coords, value = entry
			if coords[0] < b:
				if coords[1] >= b:
					ur[coords[0],coords[1]-b] = value
				else:
					ul[coords[0],coords[1]] = value
			else:
				if coords[1] >= b:
					lr[coords[0]-b,coords[1]-b] = value
				else:
					ll[coords[0]-b,coords[1]] = value
		return ur,lr,ll,ul


# construct the A matrix represeinting Poissons equation:
# - 4 * Phi(i,j) + Phi(i+1,j) + Phi(i-1,j) + Phi(i,j+1) + Phi(i,j-1) = C * rho(i,j)
def build_A(N):
	#initialize as sparse matrix
	_A = dok_matrix((N*N,N*N),dtype=sp.float32)
	
	#initialize the index array
	index_array = np.arange(N*N,dtype=np.int)
	index_array = np.reshape(index_array,(N,N))

	for i in range(N):
		for j in range(N):
			center = index_array[i,j]
			#set diagonal to -4
			_A[center,center] = -4
			#loop over the 4 non-diagonal entries
			for delta in ((-1,0),(1,0),(0,1),(0,-1)):
				_i = i + delta[0]
				_j = j + delta[1]
				#periodic wrapping
				_i = _i % N
				_j = _j % N
				index = index_array[_i,_j]
				_A[center,index] = 1
	return _A

def build_A_redblack(N):
	#initialize as sparse matrix
	_A = dok_matrix((N*N,N*N),dtype=sp.float32)
	
	#initialize the index array
	index_array = np.zeros((N,N),dtype=np.int)
	
	x = np.arange(N)
	y = np.arange(N)
	xx,yy = np.meshgrid(x,y)
	
	
	red_len = index_array[xx%2 == yy%2].size
	black_len = index_array[xx%2 != yy%2].size
	index_array[xx%2 == yy%2] = np.arange(0,red_len)
	index_array[xx%2 != yy%2] = np.arange(0,black_len) + red_len

	for i in range(N):
		for j in range(N):
			center = index_array[i,j]
			#set diagonal to -4
			_A[center,center] = -4
			#loop over the 4 non-diagonal entries
			for delta in ((-1,0),(1,0),(0,1),(0,-1)):
				_i = i + delta[0]
				_j = j + delta[1]
				#periodic wrapping
				if _i != _i%N or _j != _j % N:
					continue

				index = index_array[_i,_j]
				_A[center,index] = 1
	
	return _A,index_array


# This function improves the solution Ax=b with one Jacobi step and returns the result.
def jacobi_step(x,b,N):
	_x = x.flatten()
	_b = b.flatten()
	xnew = D.dot(_b) + DLU.dot(_x)
	return np.reshape(xnew,(N,N))

# This function improves the solution Ax=b with one Gauss-Seidel step and returns the result.
def gauss_step(x,b,N):
	_x = x.flatten()
	_b = b.flatten()
	xnew = np.array(_x,copy=True)
	offset = DU.dot(_x) + D.dot(_b)
	for i in range(_x.size):
		xnew[i] = DL.getrow(i)*_x + offset[i]
	return np.reshape(xnew,(N,N))

def gauss_step_rb(x,_b,N):
	_x = x.flatten()
	
	s = N*N/2 #scale
	
	#first do red cells
	_x[:s] =  D_r_inv.dot( -U_r.dot(_x[s:]) + _b[:s])
	_x[s:] =  D_b_inv.dot( -L_b.dot(_x[:s]) + _b[s:])

	return _x.reshape((N,N))


# This function calculates the resdiuum vector res = b - Ax, for input vectors
# of length N. The output is stored in res.
def calc_residuum(x, b, N):
	_x = x.flatten()
	_b = b.flatten()
	res = _b - A.dot(_x)
	res = np.reshape(res,(N,N))
	return res

def calc_residuum_rb(x, b, N):
	_x = x.flatten()
	res = b - A.dot(_x)
	res = np.reshape(res,(N,N))
	return res


# This function calculates the norm of the vector of length N, 
# defined as the usual quadratic vector norm.   
def norm_of_residual(res, N):
    sum = np.sum(res**2)
    return np.sqrt(sum)
    
def iteration_scheme(N):
	#defining constants
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
		#phi = jacobi_step(phi,b,N)
		phi = gauss_step(phi,b,N)

		res = calc_residuum(phi, b, N)

		r = norm_of_residual(res, N)

		print("iter=%d:  residual=%g\n"%(i, r))
		f.write("%g\n"%r)
		residuals[i] = r

	f.close()

	plt.semilogy(residuals, label="jacobi")
	plt.show()

def iteration_scheme_rb(N,index_array):
	#defining constants
	steps = 200
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
	
	b = np.arange(N*N).reshape((N,N))
	i = index_array.flatten()
	b_ = b.flatten()[i]


	for i in np.arange(steps):
		#phi = jacobi_step(phi,b,N)
		phi = gauss_step_rb(phi,b_,N)
		res = calc_residuum_rb(phi,b_, N)

		r = norm_of_residual(res, N)

		print("iter=%d:  residual=%g\n"%(i, r))
		f.write("%g\n"%r)
		residuals[i] = r

	f.close()

	plt.semilogy(residuals, label="jacobi")
	plt.show()


redblack = True

if redblack:
	N = 32

	A,index_array = build_A_redblack(N)
	U_r, D_b, L_b, D_r = to_quadrants(A)
	
	D_r_inv = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
	D_r_inv.setdiag( 1/D_r.diagonal() )
	
	D_b_inv = dok_matrix((N*N/2,N*N/2),dtype=sp.float32)
	D_b_inv.setdiag( 1/D_b.diagonal() )

	iteration_scheme_rb(N,index_array)

else:
	N = 32
	#creating the matrices as sparse matrices
	#create the Matrice A representing Poissons equation
	A = build_A(N)
	#create the inverse diagonal matrix
	d = 1./A.diagonal()
	D = dok_matrix((N*N,N*N),dtype=sp.float32)
	D.setdiag(d)
	#create the negative lower triagonal matrix
	L = -1.*sp.sparse.tril(A,-1)
	#create the negative upper triagonal matrix
	U = -1.*sp.sparse.triu(A,1)
	#create L + U
	L_U = L + U
	#create D(L + U)
	DLU = D*L_U
	#create D*U
	DU = D*U
	#create D*L as sparse row matrix
	DL = (D*L).tocsr()
	#executing the iteration scheme
	iteration_scheme(N)
