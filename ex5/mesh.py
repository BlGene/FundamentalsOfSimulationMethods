import numpy as np
def buildA(N):
    A = np.zeros((N*N,N*N))
    
    index_array = np.arange(N*N)
    index_array = index_array.reshape((N,N))
    

    # .---->y
    # |
    # |
    # v
    # x
    for x in range(N):
        for y in range(N):
            center = index_array[x,y]
            A[center,center] = -4
            
            #iterate over cell that is top,right,bottom,left
            #of center
            for delta in ((-1,0),(0,1),(1,0),(0,-1)):
                x_ = x + delta[0]
                y_ = y + delta[1]

                x_ = x_ % N
                y_ = y_ % N

                var = index_array[x_,y_]
                
                #variation right
                A[center,var] = 1
    return A




if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan)

    N = 256
    A = buildA(N)
    A

    D_array = 1/np.diag(A)
    D_matrix = np.diag(D_array)
    
    L = np.tril(A,-1)
    U = np.triu(A,1)

    x = np.arange(N*N)
    b = np.arange(N*N)
    
    for _ in range(30):
        x = np.dot(D_array,b)  + np.dot( np.dot(D_matrix,L+U), x)

    import pdb
    pdb.set_trace()

