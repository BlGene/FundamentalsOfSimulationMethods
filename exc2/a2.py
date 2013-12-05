from math import pi,sin,cos
import numpy as np
import matplotlib
import pylab as plot
import copy
import plot_frame

# defining the constants
g = 1.0
m_1 = 0.5
m_2 = 1.0
l_1 = 2.0
l_2 = 1.0

#initial values
Y_0 = [(50.)/180*pi,(-120.)/180*pi,0,0] #phi_1, phi_2, q_1, q_2
#timestep
dt = 0.05
#total number of timesteps
T_total = 2000 

#function to calculate the total energy of the system E = T + V from ph1_1, phi_2, phi_1_dot, phi_2_dot
def energy(phi_1,phi_2,phi_1_dot,phi_2_dot):
	#kinetic energy T
	T = 0.5*m_1*(l_1*phi_1_dot)**2 + 0.5*m_2*((l_1*phi_1_dot)**2+(l_2*phi_2)**2+2*l_1*l_2*phi_1_dot*phi_2_dot*cos(phi_1-phi_2))
	#potential energy V
	V = m_1*g*l_1*(1-cos(phi_1))+m_2*g*(l_1*(1-cos(phi_1))+l_2*(1-cos(phi_2)))
	return T + V

#the 2nd order differential equations in phi, phi_dot are transformed to a vectorial 1st order differential equation in Y = (phi_1,phi_2,q_1,q_2)  
#function returnig dY/dt = f(Y)
def func(Y):
	# initializing the returnvalue
	dY_ret = np.array([0.,0.,0.,0.])
	# calculating phi_1_dot
	dY_ret[0] = (l_2 * Y[2] - l_1 * cos( Y[0] - Y[1]) * Y[3])/((l_1*l_2)**2*(m_1+m_2*sin(Y[0]-Y[1])**2))
	# calculating phi_2_dot
	dY_ret[1] = ( Y[3]/m_2 - l_1 *l_2*cos( Y[0] - Y[1]) * dY_ret[0])/l_2**2
	# calculating q_1_dot
	dY_ret[2] = -1*m_2*l_1*l_2*dY_ret[0]*dY_ret[1]*sin(Y[0] - Y[1]) - (m_1 + m_2)*g*l_1*sin(Y[0])
	# calculating q_2_dot
	dY_ret[3] =    m_2*l_1*l_2*dY_ret[0]*dY_ret[1]*sin(Y[0] - Y[1]) - m_2*g*l_2*sin(Y[1])
	return dY_ret

#single step of the 2nd-order runge-kutta scheme
def step_2(Y,dt):
	k_1 = func(Y)
	k_2 = func(Y + k_1*dt)
	Y += 0.5*(k_1 + k_2)*dt
	#calculationg phi_dot to calcutate energys later
	phi_1_dot = l_2 * Y[2] - l_1 * cos( Y[0] - Y[1]) * Y[3]
	phi_2_dot = ( Y[3]/m_2 - l_1 *l_2*cos( Y[0] - Y[1]) * phi_1_dot)/l_2**2
	return Y, phi_1_dot, phi_2_dot

#single step of the 4th-order runge-kutta scheme
def step_4(Y,dt):
	k_1 = func(Y)
	k_2 = func(Y+k_1*0.5*dt)
	k_3 = func(Y+k_2*0.5*dt)
	k_4 = func(Y+k_3*dt)
	Y += (k_1/6+k_2/3+k_3/3+k_4/6)*dt
	#calculationg phi_dot to calcutate energys later
	phi_1_dot = l_2 * Y[2] - l_1 * cos( Y[0] - Y[1]) * Y[3]
	phi_2_dot = ( Y[3]/m_2 - l_1 *l_2*cos( Y[0] - Y[1]) * phi_1_dot)/l_2**2
	return Y, phi_1_dot, phi_2_dot

#complete runge-kutta-scheme. by specifying "step" 2nd- or 4th- order is chosen 
#the record flag enables saving figures to generate a movie later
def rk(step,record):
	#initializing the array for phi and q, the time and calculating the starting energy
	Y = np.array(Y_0)
	t = 0
	E_0 = energy(Y[0],Y[1],Y[2],Y[3])
	# lists to save the results of Y and E and the time
	Y_res = []
	t_res = []
	E_res = []
	Y_res.append(Y)
	t_res.append(t)
	E_res.append(E_0)
	#performing the runge-kutta-scheme for T-total times
	for i in range(T_total):
		Y,phi_1_dot,phi_2_dot = step(Y,dt)
		t += dt
		#saving results of the current step
		Y_res.append( copy.deepcopy(Y))
		t_res.append(t)
		E_res.append(energy(Y[0],Y[1],phi_1_dot,phi_2_dot))
		#if record is true save a picture of the current step
		if(record == True):
			temp = np.array(Y_res)
			plot_frame.frame(Y[0],Y[1],temp[:,0],temp[:,1],t,i)
	return Y_res,t_res,E_res

#part c: runge kutta scheme 2nd order, plotting the energy error and the angles
def part_c():
	#executing the 2nd order runge-kutta-scheme
	Y_res,t_res,E_res = rk(step_2,False)
	#casting the list to an array	
	res = np.array(Y_res)
	#initializing the figure	
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	#plotting the angles
	plot.plot(t_res,res[:,0])
	plot.plot(t_res,res[:,1])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Phi')
	plot.savefig("2nd_order_angles")
	plot.close()
	#plotting the momenta
	plot.plot(t_res,res[:,2])
	plot.plot(t_res,res[:,3])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('q')
	plot.savefig("2nd_order_momenta")
	plot.close()
	#plotting the energy error
	res = np.array(E_res)
	plot.plot(t_res,(res-res[0])/res[0])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Energy-Error')
	plot.savefig("2nd_order_deltaE")
	plot.close()

#part d: runge kutta scheme 4th order, plotting the energy error, the angle s and the generalized momenta
def part_d():
	Y_res,t_res,E_res = rk(step_4,False)
	#casting the list to an array	
	res = np.array(Y_res)
	#initializing the figure	
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	#plotting the angles
	plot.plot(t_res,res[:,0])
	plot.plot(t_res,res[:,1])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Phi')
	plot.savefig("4th_order_angles")
	plot.close()
	#plotting the momenta
	plot.plot(t_res,res[:,2])
	plot.plot(t_res,res[:,3])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('q')
	plot.savefig("4th_order_momenta")
	plot.close()
	#plotting the energy error
	res = np.array(E_res)
	plot.plot(t_res,abs(res-res[0])/res[0])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Energy-Error')
	plot.savefig("4th_order_deltaE")
	plot.close()

#part e: preparing the frames for rendering the movie via ffmpeg
def part_e():
	rk(step_2,True)


part_c()
#part_d()
#part_e()
