from math import pi,sin,cos
import numpy as np
import matplotlib
import pylab as plot
import copy

Y_0 = [(50.)/180*pi,(-120.)/180*pi,0,0] #phi_1, phi_2, q_1, q_2
dt = 0.05
T_total = 2000 

g = 1.0
m_1 = 0.5
m_2 = 1.0
l_1 = 2.0
l_2 = 1.0

def energy(phi_1,phi_2,punkt_1,punkt_2):
	T = 0.5*m_1*(l_1*punkt_1)**2 + 0.5*m_2*((l_1*punkt_1)**2+(l_2*phi_2)**2+2*l_1*l_2*punkt_1*punkt_2*cos(phi_1-phi_2))
	V = m_1*g*l_1*(1-cos(phi_1))+m_2*g*(l_1*(1-cos(phi_1))+l_2*(1-cos(phi_2)))
	return T + V

def func(Y,dt):
	dY_ret = np.array([0.,0.,0.,0.])

	dY_ret[0] = (l_2 * Y[2] - l_1 * cos( Y[0] - Y[1]) * Y[3])/((l_1*l_2)**2*(m_1+m_2*sin(Y[0]-Y[1])**2))
	dY_ret[1] = ( Y[3]/m_2 - l_1 *l_2*cos( Y[0] - Y[1]) * dY_ret[0])/l_2**2
	dY_ret[2] = -1*m_2*l_1*l_2*dY_ret[0]*dY_ret[1]*sin(Y[0] - Y[1]) - (m_1 + m_2)*g*l_1*sin(Y[0])
	dY_ret[3] =    m_2*l_1*l_2*dY_ret[0]*dY_ret[1]*sin(Y[0] - Y[1]) - m_2*g*l_2*sin(Y[1])
	return dY_ret


def step2(Y,dt):
	k_1 = func(Y,dt)
	k_2 = func(Y + k_1*dt,dt)
	Y += 0.5*(k_1 + k_2)*dt 
	punkt_1 = l_2 * Y[2] - l_1 * cos( Y[0] - Y[1]) * Y[3]
	punkt_2 = ( Y[3]/m_2 - l_1 *l_2*cos( Y[0] - Y[1]) * punkt_1)/l_2**2
	return Y, punkt_1, punkt_2

def rk_2():
	Y = np.array(Y_0)
	t = 0
	E_0 = energy(Y[0],Y[1],Y[2],Y[3])
	# lists to save the results of the temperature and the corresponding time
	Y_res = []
	t_res = []
	E_res = []
	Y_res.append(Y)
	t_res.append(t)
	E_res.append(energy(Y[0],Y[1],Y[2],Y[3]))
	i = 0
	for _ in range(T_total):
		Y,punkt_1,punkt_2 = step2(Y,dt)
		t += dt
		Y_res.append( copy.deepcopy(Y))
		t_res.append(t)
		if i < 10:
			print Y
		E_res.append(energy(Y[0],Y[1],punkt_1,punkt_2))
		i += 1
	
	return Y_res,t_res,E_res

def step_4(Y,dt):

def part_c():
	Y_res,t_res,E_res = rk_2()
	
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	
	res = np.array(Y_res)

	'''
	plot.plot(t_res,res[:,0])
	plot.plot(t_res,res[:,1])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Phi')
	plot.show()
	
	'''
	res = np.array(E_res)
	plot.plot(t_res,abs(res-res[0])/res[0])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('time')
	ax1.set_ylabel('Energy-Error')
	plot.show()

part_c()
