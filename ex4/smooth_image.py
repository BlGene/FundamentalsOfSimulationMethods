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
