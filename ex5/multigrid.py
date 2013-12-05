import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.sparse import dok_matrix

# construct the A matrix represeinting Poissons equation:
# - 4 * Phi(i,j) + Phi(i+1,j) + Phi(i-1,j) + Phi(i,j+1) + Phi(i,j-1) = C * rho(i,j)
def build_A(N):
	#initialize as sparse matrix
	_A = dok_matrix((N*N,N*N),dtype=sp.float32)
	
	#initialize the index array
	index_array = np.arange(N*N)
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
	fine_extended = np.zeros((N+1,N+1))
	return _A

# This function improves the solution Ax=b with one Gauss-Seifel step, with the
# result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
def gaussseidel_step(x, b, N, D, DU, DL):
	_x = x.flatten()
	_b = b.flatten()
	xnew = np.array(_x,copy=True)
	offset = DU.dot(_x) + D.dot(_b)
	for i in range(_x.size):
		xnew[i] = DL.getrow(i)*_x + offset[i]
	return np.reshape(xnew,(N,N))

# This function calculates the resdiuum vector res = b - Ax, for input vectors
# of length N. The output is stored in res.
def calc_residuum(x, b, N, A):
	_x = x.flatten()
	_b = b.flatten()
	res = _b - A.dot(_x)
	res = np.reshape(res,(N,N))
	return res

# This function calculates the norm of the vector of length N,
# defined as the usual quadratic vector norm.
def norm_of_residual(res, N):
    sum = np.sum(res**2)
    return np.sqrt(sum)


#This function restricts the NxN mesh stored in 'fine[ ]' to the NNxNN mesh stored in 'coarse[ ]'
def do_restrict(N, fine, NN):
	kernel = np.array(([1/16.,1/8.,1/16.],[1/8.,1/4.,1/8.],[1/16.,1/8.,1/16.]))
	fine = signal.convolve2d(fine, kernel, mode='same', boundary = 'wrap')
	coarse = fine[::2,::2]
	return coarse

# This function interpolates the the NNxNN mesh stored in 'coarse[ ]' to NxN mesh stored in 'fine[ ]'
def do_prolong(NN, coarse, N):
	coarse_pro = np.zeros((N,N))
	coarse_pro[::2,::2] = coarse[::,::]
	kernel = np.array(([1/4.,1/2.,1/4.],[1/2.,1.,1/2.],[1/4.,1/2.,1/4.]))
	fine = signal.convolve2d(coarse_pro, kernel, mode='same', boundary = 'wrap')
	return fine 

# This function carries out a V-cycle using the multigrid approach, down to N=4.
def do_v_cycle(x, b, N):
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
	x = gaussseidel_step(x, b, N, D, DU, DL)
	if N > 4:
		NN = N / 2;
		    
		#allocate some storage for the residual and error on the current and a coarsened mesh */
		err = np.zeros((N,N))
		err_coarse = np.zeros((NN, NN))
		
		#calculate the residuum
		res = calc_residuum(x, b, N, A)
		
		#restrict the residuum
		res_coarse = do_restrict(N, res, NN)
		    
		#now multiply the residuum on the coarser mesh (our b)
		#with a factor of 4 to take care of the (h) -> (2h) change, and the fact that we do not change A
		res_coarse[:,:] *= 4
		    
		#call the V-cycle again, recursively
		err_coarse = do_v_cycle(err_coarse, res_coarse, NN)
		    
		#now interpolate the error to the finer mesh
		err = do_prolong(NN, err_coarse, N)
		    
		#finally add the error to the current solution
		x[:,:] += err
    
	x = gaussseidel_step(x, b, N, D, DU, DL)
	return x

def multigrid_scheme():
	N = 256
	steps = 2000
	L = 1.0
	h = L / N
	eta = 0.1 * L
	rho0 = 10.0
	  
	#creating the matrices as sparse matrices
	#create the Matrice A representing Poissons equation
	A = build_A(N)
	
	res = np.zeros((N,N))
	
	#now set-up the density field
	x = (np.arange(-N/2,N/2)+0.5)*h
	mesh = np.meshgrid(x,x)
	r2 = mesh[0] * mesh[0] + mesh[1] * mesh[1]
	
	rho = rho0* np.exp(-r2/(2*eta**2))
	sum = np.sum(rho)
	
	#initialize the starting values for phi[] and b[]
	rho -= sum/(N**2)
	b = 4 * np.pi*h**2*rho
	phi = np.zeros((N,N))
	
	#open a file for outputting the residuum values, and then do 2000 Jacobi steps
	residuals = np.zeros(steps)
	f = file("res_multigrid.txt", "w")
	f.write("%d\n"%steps)
	
	for i in np.arange(steps):
	    phi = do_v_cycle(phi, b, N)
	    
	    res = calc_residuum(phi, b, N, A)
	    
	    r = norm_of_residual(res, N)
	    
	    print("iter=%d:  residual=%g\n"%(i, r))
	    f.write("%g\n"%r)
	    residuals[i] = r
	
	f.close()
	
	plt.semilogy(residuals, label="multigrid")
	plt.show()

multigrid_scheme()

#N = 8
#NN = N/2
#
#_coarse = np.ones((NN,NN))
#
##print _fine
##_coarse = do_restrict(N,_fine,NN)
##print _coarse
#
#print _coarse
#_fine = do_prolong(NN,_coarse,N)
#print _fine
