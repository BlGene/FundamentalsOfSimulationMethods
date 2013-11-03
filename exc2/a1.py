import matplotlib
import pylab as plot

#definition of constants
k_b = 1.38*10**-16
Lamda_0 = 10**-22 
T_0 = 20000
alpha = 10.0
beta = -0.5
n_H = 1.0
dt = 10**10
T_init = 10**7
T_end = 6000
eps_0 = 10.

#evaluate dT/dt as function of T
def func(T):
	rT = T / T_0
	if rT <= 1:
		Lamda = Lamda_0 * rT ** alpha
	else:
		Lamda = Lamda_0 * rT ** beta
	return -2./(3 * k_b ) * n_H * Lamda

#one step of the 2nd order RK-scheme
def step(T,dt):
	k_1 = func(T)
	k_2 = func(T + k_1*dt)
	return T + 0.5*(k_1 + k_2)*dt
	
#executing the RK-scheme with a fixed timestep, returning lists with the Temeperatures and corresponding times
def rk_fixed(dt):
	T = T_init
	t = 0
	# lists to save the results of the temperature and the corresponding time
	T_res = []
	t_res = []
	T_res.append(T)
	t_res.append(t)
	while T > T_end:
		T = step(T,dt)
		t += dt
		T_res.append(T)
		t_res.append(t)
	return T_res, t_res

#executing the RK-scheme with an adaptive timestep
def rk_adapt(dt_init):
	T = T_init
	t = 0
	dt = dt_init
	# list to save the results of the temperature and the corresponding time
	T_res = []
	t_res = []
	T_res.append(T)
	t_res.append(t)
	while T > T_end:
		T_a = step(T,dt)
		T_aux = step(T,dt/2.)
		T_b = step(T_aux,dt/2.)
		eps = abs(T_a - T_b)
		if eps > eps_0:
			while(eps > eps_0):
				dt /= 2.
				T_a = step(T,dt)
				T_aux = step(T,dt/2.)
				T_b = step(T_aux,dt/2.)
				eps = abs(T_a - T_b)
			T = T_b
			T_res.append(T)
			t += dt
			t_res.append(t)
		elif eps*(2**3) > eps_0:
			dt *= 2.
			T = T_b
			T_res.append(T)
			t += dt
			t_res.append(t)
		elif eps < eps_0:
			T = T_b
			T_res.append(T)
			t += dt
			t_res.append(t)
	
	return T_res,t_res


# RK-scheme second order with fixed timestep dt = 10**10 
# Plot with Logarithmic t axis and linear T axis
def part_a():
	dt = 10**10
	T_res,t_res = rk_fixed(dt)
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	plot.plot(t_res,T_res)
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xscale('log')
	ax1.set_xlabel('time 10**10sec')
	ax1.set_ylabel('Temp K')
	plot.show()

# RK-scheme second order with fixed timestep for different timesteps dt
def part_b():
	j = 0
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	for i in (0.8,1.0,2.0,5.0,10.,15,20):
		dt = i * 10**10
		T_res,t_res = rk_fixed(dt)
		
		plot.plot(t_res,T_res,label=str(i))
		handles, labels = ax1.get_legend_handles_labels()
		ax1.legend(handles, labels)
		ax1.set_xscale('log')
		ax1.set_xlabel('time 10**10sec')
		ax1.set_ylabel('Temp K')
		
		j += 1
	plot.show()

def part_c():
	j = 0
	fig1 = plot.figure()
	ax1 = fig1.add_subplot(111)
	for i in (0.8,1.0,10,20.,30.0,100,1000):
		dt = i * 10**10
		T_res,t_res = rk_adapt(dt)
		print len(T_res),len(t_res)	
		plot.plot(t_res,T_res,label=str(i))
		handles, labels = ax1.get_legend_handles_labels()
		ax1.legend(handles, labels)
		ax1.set_xscale('log')
		ax1.set_xlabel('time 10**10sec')
		ax1.set_ylabel('Temp K')
		
		j += 1
	plot.show()
	
part_c()
