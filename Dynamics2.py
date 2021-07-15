import numpy as np
import math as ma
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft,ifft
import time
# find peak value of displacement, veloctiy, total acceleration for different periods.  
start = time.process_time()

filename = 'elcentro.dat'
indata = np.loadtxt(filename, usecols=(0,1))
time = (indata[:,0])
excit_F = np.multiply((indata[:,1]), 9.8)
# earth_data = np.sin(time)
earth_data = np.multiply((indata[:,1]), 9.8)

# define paremeter

def trange(start, final, dt):
	trange = np.arange(start, final, dt)
	# print(len(trange), 'with', dt,'and', final)
	return trange

def compute_amp(u_c, u_cd, u_cdd, dt, finaltime, c, k, beta, gamma):
	timerange = np.arange(0, finaltime, dt)
	amp = [u_c]
	velocity = [u_cd]
	acceleration = [u_cdd]
	for i in range(0,len(timerange)-1):
		u0_np1 = u_c + dt*u_cd + 0.5*(dt**2)*(1-2*beta)*u_cdd
		u0_dnp1 = u_cd + dt*(1 - gamma)*u_cdd

		u_ddnp1	=(excit_F[i] - c * u0_dnp1 - k*u0_np1)/(1 + gamma * c*dt + beta * k * np.square(dt))
		
		u_np1 = u_c + dt*u_cd + 0.5*(dt**2)*((1-2 * beta) * u_cdd + 2*beta*u_ddnp1)
		u_dnp1 = u_cd + dt*((1 - gamma)* u_cdd + gamma*u_ddnp1)

		amp.append(u_np1)
		velocity.append(u_dnp1)
		acceleration.append(u_ddnp1)
		u_c, u_cd, u_cdd = u_np1,u_dnp1, u_ddnp1
	
	return amp,	velocity, acceleration
def abc(T):
	
	
	# for k in r:
	frequency = (2*np.pi)/T
	beta = 1/4
	gamma = 1/2
	dt = 0.02
	m = 1000
	omega = frequency
	chi = 0.05
	k = np.square(omega)
	c = chi * omega
	# dk = k + (gamma/(beta*dt)) * c + (1/ (beta* np.square(dt)))* m
	u_0 = 0
	ud_0 = 0
	udd_0 = 0
	
	displacement, veloctiy, acceleration = compute_amp(u_0, ud_0, udd_0, dt, 31.18, c, k, beta, gamma)
	
	
	dis = max(displacement)
	vel = max(veloctiy)
	acc =max(acceleration)
	return dis, vel, acc, displacement, veloctiy, acceleration

# r = np.arange(10, 0.1 ,-0.1)
max_dis= []
max_vel = []
max_acc = []
analytical_max_vel = []
analytical_max_acc = []
T = np.arange(0.01,10,0.05)
for DT in T:
	frequency = (2*np.pi)/DT
	dis = abc(DT)[0]
	vel = abc(DT)[1]
	acc = abc(DT)[2]
	max_dis.append(dis)
	max_vel.append(vel)
	max_acc.append(acc)
	analytical_max_vel.append(dis*frequency)
	analytical_max_acc.append(dis*frequency*frequency)
np.save('period.npy', T)
np.save('acc_spectrum.npy', max_acc)
np.save('dis_spectrum.npy', max_dis)
displacement = abc(0.5)[3]
velocity = abc(0.5)[4]
acceleration = abc(0.5)[5]
chi = 0.05


plt.figure(1)
plt.plot(max_dis)
plt.grid()
plt.xlabel('period steps(s)')

plt.ylabel('max Displacement')
plt.title('displacement spectra - '+ r'$\xi = $'+str(chi))
plt.savefig('displacement spectra')

plt.figure(2)
plt.plot(max_vel,label= 'pseudo velocity spectra')
plt.plot(analytical_max_vel, label='approximate pseudo velocity spectra')
plt.xlabel('period steps(s)')

plt.ylabel('max velocity')
plt.title('velocity spectra - '+ r'$\xi = $'+str(chi))
plt.grid()
plt.legend()
plt.savefig('velocity spectra')

plt.figure(3)
plt.plot(max_acc,label= 'pseudo acceleration spectra')
plt.plot(analytical_max_acc, label='approximate pseudo acceleration spectra')
plt.xlabel('period steps(s)')

plt.ylabel('max acceleration')
plt.title('acceleration spectra - '+ r'$\xi = $'+str(chi))
plt.grid()
plt.legend()
plt.savefig('acceleration spectra')


# plt.figure(4)
# plt.plot(earth_data)
# plt.xlabel('time step')
# plt.ylabel('earthquake acceleration')
# plt.title('earthquake data')
# plt.grid()
# plt.savefig('earthquake data')

#######
# plt.figure(5, figsize = (15,5))

# plt.subplot(1,3,1)
# plt.plot(abc(0.5)[3], label = 'period = 0.5')
# plt.plot(abc(5)[3], label = 'period = 5')

# plt.xlabel('time step')
# plt.ylabel('displacement')
# plt.legend()
# plt.grid()

# plt.subplot(1,3,2)
# plt.plot(abc(0.5)[3], label = 'period = 0.5')

# plt.plot(abc(10)[3], label = 'period = 10')
# plt.xlabel('time step')
# plt.legend()
# plt.grid()

# plt.subplot(1,3,3)

# plt.plot(abc(5)[3], label = 'period = 5')
# plt.plot(abc(10)[3], label = 'period = 10')
# plt.xlabel('time step')
# plt.legend()
# plt.grid()

# plt.suptitle('Displacement-Time with various period'+' when '+r'$\xi = $'+str(chi))
# ######

# plt.figure(6)
# plt.plot(acceleration)
# plt.xlabel('time step')
# plt.ylabel('acceleration')
# plt.title('acceleration-time(period = 5)')
# plt.grid()

plt.show()




