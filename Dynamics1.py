import numpy as np
import math as ma
import matplotlib 
import matplotlib.pyplot as plt
from array import *
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.fftpack import fft,ifft

# coff_omega = np.arange(0.5, 10, 0.5)
# omega_f = coff_omega*ma.pi
# k = np.square(omega_f)
# m = 1
# omega = ma.sqrt(k/m)
# # print('omega =', omega)
# # dt = 0.01
# period = 2*ma.pi/omega
# # print('period =',period)
# # timerange = np.arange(0, period, dt)
u0 = 0
u0d = 2
u0dd = 0
# print('dt =', dt)
#Newmark method parameter
gama = 0.5	
#Numerical Solution
u_c, u_cd, u_cdd = u0, u0d, u0dd
ampNum1, ampNum2, ampNum3, ampNum4 = [], [], [], []
ampNum1.append(u_c)
ampNum2.append(u_c)
beta1 = 1/2
beta2 = 1/6
omega_f = 2*ma.pi
k = np.square(omega_f)
m = 1
omega = ma.sqrt(k/m)
period = 2*ma.pi/omega
# stable condition

# print('limit =' ,limit)
#
def trange(start, final, dt):
	trange = np.arange(start, final, dt)
	# print(len(trange), 'with', dt,'and', final)
	return trange

def compute_amp(beta, u_c, u_cd, u_cdd, dt, finaltime ):
	timerange = np.arange(0, finaltime, dt)
	amp = [u_c]
	velocity = [u_cd]
	for i in range(0,len(timerange)-1):
		u0_np1 = u_c + dt*u_cd + 0.5*(dt**2)*(1-2*beta)*u_cdd
		u0_dnp1 = u_cd + dt*(1 - gama)*u_cdd

		u_ddnp1	=(-k*u0_np1)/(1 + beta * k * dt**2)

		u_np1 = u_c + dt*u_cd + 0.5*(dt**2)*((1-2 * beta) * u_cdd + 2*beta*u_ddnp1)
		u_dnp1 = u_cd + dt*((1 - gama)* u_cdd + gama*u_ddnp1)

		amp.append(u_np1)
		velocity.append(u_dnp1)

		u_c, u_cd, u_cdd = u_np1,u_dnp1, u_ddnp1
	return amp,	velocity

# plotting
def amplitudeplot(start, cofficient_period, dt, beta, u_c, u_cd, u_cdd, plotorder):
	plt.subplot(3, 2, plotorder)
	plt.plot(trange(start, cofficient_period*period, dt), compute_amp(beta, u_c, u_cd, u_cdd, dt, cofficient_period*period)[0], 'b',label='dt ='+str(dt))
	plt.grid()
	# plt.title('Numerical Solutiuon with Beta = 0.5')
	plt.xlabel('time (s)')

def square(num):
	return num * num

def velo(start, cofficient_period, dt, beta, u_c, u_cd, u_cdd, plotorder):
	plt.subplot(4, 1, plotorder)
	
	plt.plot(trange(start, cofficient_period*period, dt), 0.5*(np.square(compute_amp(beta, u_c, u_cd, u_cdd, dt, cofficient_period*period)[1])), 'b',label=r'$\Delta$t = '+str(dt))
	plt.grid()
	# if dt >= 1:
	# 	plt.ylim([0,5])
	# 	plt.xlim([0,20])
	plt.legend(loc = 'upper right')
	# plt.title('Numerical Solutiuon with Beta = 0.5')

def errorplot(start, cofficient_period, dt, beta, u_c, u_cd, u_cdd, plotorder, omega):
	plt.subplot(4, 1, plotorder)
	Wave = ( u_cd/omega)*np.sin(omega*trange(start, cofficient_period*period, dt))
	plt.plot(trange(start, cofficient_period*period, dt), ((compute_amp(beta, u_c, u_cd, u_cdd, dt, cofficient_period*period)[0])), 'b',label=r'$\Delta$t = '+str(dt))
	plt.grid()
	plt.legend(loc = 'upper right')


	# plot amplitude
	# plt.figure(1)
	# amplitudeplot(0, 3, 0.1, beta1, u0, u0d, u0dd, 1)
	# amplitudeplot(0, 3, 0.3, beta1, u0, u0d, u0dd, 3)
	# amplitudeplot(0, 3, 0.6, beta1, u0, u0d, u0dd, 5)
	# amplitudeplot(0, 300, 0.1, beta1, u0, u0d, u0dd, 2)
	# amplitudeplot(0, 300, 0.3, beta1, u0, u0d, u0dd, 4)
	# amplitudeplot(0, 300, 0.6, beta1, u0, u0d, u0dd, 6)

	# plot velocity
	# beta = 1/2
	# fig1 = plt.figure(1, figsize = (15,10))
	# velo(0, 100, 0.25, beta1, u0, u0d, u0dd, 1)
	# velo(0, 100, 0.5, beta1, u0, u0d, u0dd, 2)
	# velo(0, 100, 1, beta1, u0, u0d, u0dd, 3)
	# velo(0, 100, 2, beta1, u0, u0d, u0dd, 4)
	# plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{2}$'  + '   Kinetic Energy with various '+ r'$\Delta$' + 't', fontsize = 15)
	# fig1.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
	# fig1.text(0.06, 0.5, 'Energy', ha='center', va='center', rotation='vertical')




# beta = 1/6
# fig2 = plt.figure(2, figsize = (15,10))
# velo(0, 100, 0.25, beta2, u0, u0d, u0dd, 1)
# velo(0, 100, 0.5, beta2, u0, u0d, u0dd, 2)
# velo(0, 3, 1, beta2, u0, u0d, u0dd, 3)
# velo(0, 3, 2, beta2, u0, u0d, u0dd, 4)
# plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{6}$'  + '   Kinetic Energy with various '+ r'$\Delta$' + 't', fontsize = 15)
# fig2.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
# fig2.text(0.06, 0.5, 'Energy', ha='center', va='center', rotation='vertical')

	# fig2 = plt.figure(2, figsize = (15,10))
	# velo(0, 10, 0.25, beta2, u0, u0d, u0dd, 1)
	# velo(0, 10, 0.5, beta2, u0, u0d, u0dd, 2)
	# velo(0, 10, 1, beta2, u0, u0d, u0dd, 3)
	# velo(0, 10, 2, beta2, u0, u0d, u0dd, 4)
	# plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{6}$'  + '   Kinetic Energy with various '+ r'$\Delta$' + 't', fontsize = 15)
	# fig2.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
	# fig2.text(0.06, 0.5, 'Energy', ha='center', va='center', rotation='vertical')

	#error vis
	# fig3 = plt.figure(3, figsize = (15,10))
	# errorplot(0, 100, 0.25, beta1, u0, u0d, u0dd, 1, omega)
	# errorplot(0, 100, 0.5, beta1, u0, u0d, u0dd, 2, omega)
	# errorplot(0, 100, 1, beta1, u0, u0d, u0dd, 3, omega)
	# errorplot(0, 100, 2, beta1, u0, u0d, u0dd, 4, omega)
	# plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{2}$'  + '   Error between Analytical and Numerical Solution with various '+ r'$\Delta$' + 't', fontsize = 15)
	# fig3.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
	# fig3.text(0.06, 0.5, 'Error', ha='center', va='center', rotation='vertical')


fig4 = plt.figure(4, figsize = (15,10))
errorplot(0, 100, 0.25, beta2, u0, u0d, u0dd, 1, omega)
errorplot(0, 100, 0.5, beta2, u0, u0d, u0dd, 2, omega)
errorplot(0, 5, 1, beta2, u0, u0d, u0dd, 3, omega)
errorplot(0, 5, 2, beta2, u0, u0d, u0dd, 4, omega)
plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{6}$'  + '   Error between Analytical and Numerical Solution with various '+ r'$\Delta$' + 't', fontsize = 15)
fig4.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
fig4.text(0.06, 0.5, 'Error', ha='center', va='center', rotation='vertical')

	# fig5 = plt.figure(5,figsize = (15,10))
	# errorplot(0, 10, 0.25, beta2, u0, u0d, u0dd, 1, omega)
	# errorplot(0, 10, 0.5, beta2, u0, u0d, u0dd, 2, omega)
	# errorplot(0, 10, 1, beta2, u0, u0d, u0dd, 3, omega)
	# errorplot(0, 10, 2, beta2, u0, u0d, u0dd, 4, omega)
	# plt.suptitle(r'$\beta$' + '=' + r'$\frac{1}{6}$'  + '   Partial enlarged view on Error '+ r'$\Delta$' + 't', fontsize = 15)
	# fig5.text(0.5, 0.07, 'Time(s)', ha='center', va='center')
	# fig5.text(0.06, 0.5, 'Error', ha='center', va='center', rotation='vertical')

#ff ana

# dt = np.arange(0.1,0.2, 0.01)

# for i in dt:
omega_f = 2*ma.pi
k = np.square(omega_f)
m = 1
omega = ma.sqrt(k/m)
period = 2*ma.pi/omega

# print('period =',period)
# timerange = np.arange(0, period, dt)

# beta = [1/2, 0, 1/6, 1/12, 1/3]
# u0 = 0
# u0d = 2
# u0dd = 0
### plot period shift
# plt.figure(1)

# def temporary(b, k):
# 	N_P = []
# 	dt_array = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,0.4,0.45,0.5,0.6]
# 	for dt in dt_array:
# 		t = np.arange(0, 100*period, dt)
# 		if b < 1/4:
# 			limit = ma.sqrt(1/((1/4- b)* (k))) #define the limit, program would collapse when unstable
# 			if dt > limit:
# 				continue
# 		phi = np.arctan((omega * u0) / u0d)
# 		u_max = np.sqrt(square(u0) + square(u0d/omega))
# 		# print(b, u_c, u_cd, u_cdd, dt, 1000*period)
# 		u_dispalcement = u_max * np.sin(omega*t + phi)
# 		u_c_newmarks= compute_amp(b, u_c, u_cd, u_cdd, dt, 100*period)[0]
# 		print(max(u_c_newmarks), b, dt, max(u_dispalcement))
# 		# FFT
# 		u_dispalcement_f = fft(u_dispalcement)
# 		u_c_newmarks_f = fft(u_c_newmarks)

# 		ff_analy = abs(u_dispalcement_f)
# 		ff_newmarks = abs(u_c_newmarks_f)

# 		ff_analy=abs(u_dispalcement_f)/len(t) 
# 		ff_newmarks=abs(u_c_newmarks_f)/len(t)

# 		ff_analy = ff_analy[range(int(len(t)/2))]
# 		ff_newmarks = ff_newmarks[range(int(len(t)/2))]


# 		nat_fre = int(np.where(ff_analy == max(ff_analy))[0])
# 		# print(np.where(ff_newmarks == max(ff_newmarks))[0], dt, b)
# 		newmarks_fre = int(np.where(ff_newmarks == max(ff_newmarks))[0])
# 		rate = nat_fre/newmarks_fre
# 		New_Period = period*rate
# 		N_P.append(New_Period)
# 	N_P = np.array(N_P)
# 	dT = N_P - period
	
# 	X, Y = [], []
# 	index = 0
# 	for k in dt_array:

# 		if b < 1/4:
# 			if k > limit:
# 				continue
# 		small_x = k/period
# 		small_y = dT[index]/period
# 		# print(small_y, b, k)
# 		index += 1
# 		X.append(small_x)
# 		Y.append(small_y)
# 	plt.plot(X,Y)
	
# 	return 

# for b in beta:
# 	temporary(b, k)


# plt.xlabel(r'$\Delta t/T$')
# plt.ylabel(r'$\Delta T/T$')
# plt.title('Period error with' + r'$\gamma = \frac{1}{2}$' + ' - various '+ r'$\beta$')
# plt.legend([])
# plt.legend([r'$\beta = \frac{1}{2}$', r'$\beta = 0$', r'$\beta = \frac{1}{6}$', r'$\beta = \frac{1}{12}$', r'$\beta = \frac{1}{3}$'])
# plt.xlim([0, 0.5])
# plt.grid()
### plot period shift



# plt.grid()
# plt.figure()
# plt.plot(ff_analy, label = 'Analytical Solution')
# plt.plot(ff_newmarks, label = 'Numerical Solution')
# plt.grid()
# plt.legend()
# plt.xlabel('Frequency')
# plt.ylabel('amplitude')
# plt.title('Frequency - amplitude')
plt.show()

