import math
import numpy as np
import matplotlib.pyplot as plot

#Reads a square image in 8-bit/color PPM format from the given file. Note: No checks on valid format are done.
def readImage(filename):
    f = file(filename,"rb")
    
    f.readline()
    s = f.readline()
    f.readline()
    (pixel, pixel) = [t(s) for t,s in zip((int,int),s.split())]
    
    data = np.fromfile(f,dtype=np.uint8,count = pixel*pixel*3)
    img = data.reshape((pixel,pixel,3)).astype(np.double)
    
    f.close()
    
    return img, pixel
    

#Writes a square image in 8-bit/color PPM format.
def writeImage(filename, image):
    f = file(filename,"wb")
    
    pixel = image.shape[0]
    f.writelines("P6\n%d %d\n%d\n"%(pixel, pixel, 255))
    
    image = image.astype(np.uint8)
    
    image.tofile(f)
    
    f.close()
    
    
def smooth():

    img, pixel = readImage("aq-original.ppm")


    #Now we set up our desired smoothing kernel. We'll use complex number for
    #it even though it is real. 
    kernel_real = np.zeros((pixel,pixel),dtype=np.complex)

    hsml = 10.

    #now set the values of the kernel 
    for i in np.arange(pixel):
        for j in np.arange(pixel):
            
            r = (i*i + j*j)**.5 / hsml
            
            if abs(r) < .5:
                val = 1 - 6. * r**2 + 6.* r**3
            elif abs(r) < 1.0:
                val = 2*(1-r)**3
            else:
                val = 0

            kernel_real[i, j] =  val
    
    #Now normalize the image
    kernel_real /= np.sum(kernel_real)

    #Let's calculate the Fourier transform of the kernel
    kernel_kspace = np.fft.fft2(kernel_real)

    #we now convolve each color channel with the kernel using FFTs
    for colindex in np.arange(3):
       
        #further space allocations for image transforms
        color_real = np.zeros((pixel,pixel),dtype=np.complex)
        
        #copy input color into complex array
        color_real[:,:].real = img[:,:,colindex]
        
        sum_before = np.sum(color_real)

        #forward transform
        color_kspace = np.fft.fft2(color_real)
        
        #multiply with kernel in Fourier space
        #TODO: fill in code here
        color_kspace = np.multiply(color_kspace,kernel_kspace)

        #backward transform
        color_real = np.fft.ifft2(color_kspace)
        
        sum_after = np.sum(color_real)

        print "Before:",sum_before.real,
        print "\tAfter:", sum_after.real,
        print "\tDiff:",sum_before.real-sum_after.real,
        print "\t%dev:",(sum_before.real-sum_after.real)/sum_before.real
        
        #copy real value of complex result back into color array
        img[:,:,colindex] = color_real.real
        

    writeImage("aq-smoothed.ppm", img)

'''
Output for each channel assuming RGB ordering, deviation is super small

Before: 25182443.0 	After: 25182443.0 	Diff: -3.94880771637e-07 	%dev: -1.5680796801e-14
Before: 12287845.0 	After: 12287845.0 	Diff: 5.96046447754e-08 	%dev: 4.85069959585e-15
Before: 24795222.0 	After: 24795222.0 	Diff: -5.99771738052e-07 	%dev: -2.41890045611e-14
'''


def gauss_func(dists):
    height = 1.0
    sigma = 0.07
    return np.sum( height*np.exp( -(dists/sigma)**4), axis=2)

def pot_func(dists):
    #g_x,g_y,g_z  =  np.gradient(-4*math.pi/(dists**4+.004))
    #return g_x*g_x #funky umages
    
    mesh = np.sum( -4*math.pi/(dists**2+0.001),axis=2)
    g_x,g_y  =  np.gradient(mesh,1./256,1./256)
    return -1*(g_x*g_x + g_y*g_y)**.5
    return mesh

#used for calculating the mass field
def cloud_dist(particle_centers,L):
   
    def f(x,y,function):
        exp_dim = particle_centers.shape[0]         #256x256
        
        x_exp = np.expand_dims(x,2)                 #256x256xN_particles
        y_exp = np.expand_dims(y,2)
        
        x_tiled = np.repeat(x_exp,exp_dim,axis=2)   #fill the 3rd dimension
        y_tiled = np.repeat(y_exp,exp_dim,axis=2)
        
        x_tiled -= particle_centers[:,0]            #subtract each particle center
        y_tiled -= particle_centers[:,1]            #from one of the new dimensions
        
        x_tiled = np.where(np.fabs(x_tiled)>L/2.,x_tiled-np.rint(x_tiled/L)*L,x_tiled)
        y_tiled = np.where(np.fabs(y_tiled)>L/2.,y_tiled-np.rint(y_tiled/L)*L,y_tiled)

        dists = (x_tiled*x_tiled+y_tiled*y_tiled)**.5   #compute distance for each 256x256 element

        return function(dists)                  #compute function for each 256x256 element
        
    return f

#need to consider force directions here
def force_dist(particle_centers,L):
   
    def f(x,y,function):
        exp_dim = particle_centers.shape[0]         #256x256
        
        x_exp = np.expand_dims(x,2)                 #256x256xN_particles
        y_exp = np.expand_dims(y,2)
        
        x_tiled = np.repeat(x_exp,exp_dim,axis=2)   #fill the 3rd dimension
        y_tiled = np.repeat(y_exp,exp_dim,axis=2)
        
        x_tiled -= particle_centers[:,0]            #subtract each particle center
        y_tiled -= particle_centers[:,1]            #from one of the new dimensions
        
        x_tiled -= np.rint(x_tiled/L)*L
        y_tiled -= np.rint(y_tiled/L)*L

        dists = (x_tiled*x_tiled+y_tiled*y_tiled)**.5   #compute distance for each 256x256 element

        return function(dists)                  #compute function for each 256x256 element
        
    return f


#used for constructing the kernel
def inv_dist(x,y,L):
    x_tiled = np.where(np.fabs(x)>L/2.,x-np.rint(x/L)*L,x)
    y_tiled = np.where(np.fabs(y)>L/2.,y-np.rint(y/L)*L,y)

    dists = (x_tiled*x_tiled+y_tiled*y_tiled)**2   #compute distance for each 256x256 element
    
    return 1.0/dists


def particle_mesh_force():
    L = 1.0
    N_particles = 3
    N_grid = 256
   
    particle_centers = np.c_[
        np.random.uniform(0,L,N_particles),
        np.random.uniform(0,L,N_particles) ]
    
    #some debugging
    #particle_centers = np.r_[particle_centers,np.array(((.99,.99),))]
    #particle_centers = np.r_[particle_centers,np.array(((.01,.01),))]
    
    f = cloud_dist(particle_centers,L)
    g = force_dist(particle_centers,L)

    x=np.linspace(0.0,L,N_grid)
    y=np.linspace(0.0,L,N_grid)
    
    xx,yy = np.meshgrid(x,y)
    
    mesh = f(xx,yy,gauss_func)
    
    direct = g(xx,yy,pot_func)
    
    kernel_kspace = -4*math.pi*inv_dist(xx,yy,L)
    kernel_kspace = np.where(np.isinf(kernel_kspace),0,kernel_kspace)
    kernel_kspace = np.array(kernel_kspace,dtype=complex)

    mass_real = np.array( mesh ,dtype=np.complex)
       
    #forward transform
    mass_kspace = np.fft.fft2(mass_real)
   
    #multiply with kernel in Fourier space
    mass_kspace = np.multiply(mass_kspace,kernel_kspace)
    
    mass_real = np.fft.ifft2(mass_kspace)
    
    fig = plot.figure()
    ax1 = fig.add_subplot(1,3,1)
    ax1.imshow(mesh)
    ax2 = fig.add_subplot(1,3,2)
    ax2.imshow(mass_real.real)
    ax3 = fig.add_subplot(1,3,3)
    ax3.imshow(direct)

    plot.show()
 
 
if __name__ == "__main__":
    #smooth()
    particle_mesh_force()
