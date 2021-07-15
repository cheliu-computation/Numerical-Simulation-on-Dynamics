#dynamics 3
import numpy as np
import math as ma
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft,ifft
import time
from numpy import linalg as LA
import sympy as sp
from sympy import *
from sympy.matrices import Matrix


m1 = 400
m2 = 300
m3 = 200
diemension_column = 0.35
E = 25*1000000000 # GN/M^2
height = 3
h_cubic = np.power(height, 3)
I = np.power(diemension_column, 4)/12
k = 12*(E*I)*6/h_cubic
K_matrix = np.zeros((3,3))
K_matrix[0,0] = k*2
K_matrix[1,1] = k*2
K_matrix[2,2] = k
K_matrix[0,1] = -k
K_matrix[1,0] = -k
K_matrix[2,1] = -k
K_matrix[1,2] = -k
K_matrix[0,2] = 0
K_matrix[2,0] = 0

# mass matrix
M_matrix = np.zeros((3,3))
M_matrix[0,0] = 4
M_matrix[1,1] = 3
M_matrix[2,2] = 2
M_matrix = M_matrix*100000

# finding the natrual frequency

### --- Q1
M_inverse = np.linalg.inv(M_matrix) 
matrix = np.matmul(K_matrix, M_inverse)
eig = LA.eig(matrix)
na_f = np.sqrt(eig[0])
eig_V = eig[1]
period_Q1 = 2*ma.pi/na_f
T1, T2, T3 = period_Q1[-1], period_Q1[-2], period_Q1[-3]
V1_Q1, V2_Q1, V3_Q1 = eig_V[:,-1], eig_V[:,-2], eig_V[:,-3] # modal shape vector
# load spectrum
acc_spectrum = np.load('acc_spectrum.npy')
dis_spectrum = np.load('dis_spectrum.npy')
period = np.load('period.npy')

# interplotation
def find_peak_acc(index, spectrum, T):
	x1 = period[index]
	y1 = spectrum[index]
	x3 = period[index+1]
	y3 = spectrum[index+1]
	x2 = T
	y2 = y1 + (x2-x1) * (y3 - y1)/(x3 - x1)
	return y2

left = 0
right = len(period)-1
def Bsearch(left,right,T):
	while left < right:
		mid = (left + right) // 2
		if period[mid] < T and period[mid+1] > T:
			return(mid)
			break
		elif period[mid] < T:
			left = left+1
		elif period[mid] > T:
			right = right - 1

def cal_displacement(dis, v, M):
	d = np.ones((3,1))
	m_bar = np.dot(np.dot(v.transpose(), M), v) # modal mass
	l = np.dot(np.dot(v.transpose(), M), d)/m_bar # load participation
	U = l*dis*v # displacement
	return U

index_T1 = Bsearch(left,right,T1)
index_T2 = Bsearch(left,right,T2)
index_T3 = Bsearch(left,right,T3)

###  peak displacement for question 1
peak_dis_1 = find_peak_acc(index_T1, dis_spectrum, T1)
peak_dis_2 = find_peak_acc(index_T2, dis_spectrum, T2)
peak_dis_3 = find_peak_acc(index_T3, dis_spectrum, T3)

peak_dis_M1 = cal_displacement(peak_dis_1, V1_Q1, M_matrix)
peak_dis_M2 = cal_displacement(peak_dis_2, V2_Q1, M_matrix)
peak_dis_M3 = cal_displacement(peak_dis_3, V3_Q1, M_matrix)
print('Peak displacement vector of mode1 = ' + str(peak_dis_M1))
print('Peak displacement vector of mode2 = ' + str(peak_dis_M2))
print('Peak displacement vector of mode3 = ' + str(peak_dis_M3))
### 
peak_acc_1 = find_peak_acc(index_T1, acc_spectrum, T1)
peak_acc_2 = find_peak_acc(index_T2, acc_spectrum, T2)
peak_acc_3 = find_peak_acc(index_T3, acc_spectrum, T3)

sum_MZ = (M_matrix[0,0]*3 + M_matrix[1,1]*6 + M_matrix[2,2]*9)
sum_M = M_matrix[0,0]+M_matrix[1,1]+M_matrix[2,2]

V_b1 = sum_M*peak_acc_1
V_b2 = sum_M*peak_acc_2
V_b3 = sum_M*peak_acc_3
## base shear force for three modes
F1_Vb1, F2_Vb1, F3_Vb1 = V_b1 * (M_matrix[0,0]*3)/sum_MZ, V_b1 * (M_matrix[1,1]*6)/sum_MZ, V_b1 * (M_matrix[2,2]*9)/sum_MZ
F_M1 = [F1_Vb1, F2_Vb1, F3_Vb1]

F1_Vb2, F2_Vb2, F3_Vb2 = V_b2 * (M_matrix[0,0]*3)/sum_MZ, V_b2 * (M_matrix[1,1]*6)/sum_MZ, V_b2 * (M_matrix[2,2]*9)/sum_MZ
F_M2 = [F1_Vb2, F2_Vb2, F3_Vb2]

F1_Vb3, F2_Vb3, F3_Vb3 = V_b3 * (M_matrix[0,0]*3)/sum_MZ, V_b3 * (M_matrix[1,1]*6)/sum_MZ, V_b3 * (M_matrix[2,2]*9)/sum_MZ
F_M3 = [F1_Vb3, F2_Vb3, F3_Vb3]
print('Base shear force for Mode1 = ' + str(F_M1))
print('Base shear force for Mode2 = ' + str(F_M2))
print('Base shear force for Mode3 = ' + str(F_M3))

#SRSS
srss_dis = np.sqrt(np.square(peak_dis_M1)+np.square(peak_dis_M2)+np.square(peak_dis_M3))
srss_acc = np.sqrt(np.square(F_M1)+np.square(F_M2)+np.square(F_M3))
print('srss dis = '+str(srss_dis))
print('srss acc = '+str(srss_acc))

## base shear force for three modes

###### Question 2
# Ritz Vector
r1 = [[1], [2], [3]]
r2 = [[1], [4], [9]]
r1 = np.array(r1)
r2 = np.array(r2)
R = np.hstack((r1,r2))
## M_hat and K_hat 
M_hat = np.matmul(np.matmul(R.transpose(), M_matrix), R)
K_hat = np.matmul(np.matmul(R.transpose(), K_matrix), R)
## eig value problem
M_inverse_Q2 = np.linalg.inv(M_hat) 
matrix_Q2 = np.matmul(K_hat, M_inverse_Q2)
eig_Q2 = LA.eig(matrix_Q2)
na_f_Q2 = np.sqrt(eig_Q2[0])
eig_V_Q2 = eig_Q2[1]

period_Q2 = 2*ma.pi/na_f_Q2
T1_Q2, T2_Q2 = period_Q2[0], period_Q2[1]
X1_Q2, X2_Q2 = eig_V_Q2[:,0], eig_V_Q2[:,1]

V1_Q2, V2_Q2 = np.dot(R,X1_Q2), np.dot(R,X2_Q2)

index_T1_Q2 = Bsearch(left,right,T1_Q2)
index_T2_Q2 = Bsearch(left,right,T2_Q2)

###  peak displacement for question 2
def cal_displacement_Q2(dis, v, M):
	d = np.ones((3,1))
	m_bar = np.dot(np.dot(v.transpose(), M), v) # modal mass
	l = np.dot(np.dot(v.transpose(), M), d)/m_bar # load participation
	U = l*dis*v # displacement
	return U

peak_dis_1_Q2 = find_peak_acc(index_T1_Q2, dis_spectrum, T1_Q2)
peak_dis_2_Q2 = find_peak_acc(index_T2_Q2, dis_spectrum, T2_Q2)

peak_dis_M1_Q2 = cal_displacement_Q2(peak_dis_1_Q2, V1_Q2, M_matrix)
peak_dis_M2_Q2 = cal_displacement_Q2(peak_dis_1_Q2, V2_Q2, M_matrix)
print('						')
print('Peak displacement vector of mode1(Q2) = ' + str(peak_dis_M1_Q2))
print('Peak displacement vector of mode2(Q2) = ' + str(peak_dis_M2_Q2))

peak_acc_1_Q2 = find_peak_acc(index_T1_Q2, acc_spectrum, T1_Q2)
peak_acc_2_Q2 = find_peak_acc(index_T2_Q2, acc_spectrum, T2_Q2)

sum_MZ = (M_matrix[0,0]*3 + M_matrix[1,1]*6 + M_matrix[2,2]*9)
sum_M = M_matrix[0,0]+M_matrix[1,1]+M_matrix[2,2]

V_b1_Q2 = sum_M*peak_acc_1_Q2
V_b2_Q2 = sum_M*peak_acc_2_Q2

## base shear force for three modes
F1_Vb1_Q2, F2_Vb1_Q2, F3_Vb1_Q2 = V_b1_Q2 * (M_matrix[0,0]*3)/sum_MZ, V_b1_Q2 * (M_matrix[1,1]*6)/sum_MZ, V_b1_Q2 * (M_matrix[2,2]*9)/sum_MZ
F_M1_Q2 = [F1_Vb1_Q2, F2_Vb1_Q2, F3_Vb1_Q2 ]

F1_Vb2_Q2, F2_Vb2_Q2, F3_Vb2_Q2 = V_b2_Q2 * (M_matrix[0,0]*3)/sum_MZ, V_b2_Q2 * (M_matrix[1,1]*6)/sum_MZ, V_b2_Q2 * (M_matrix[2,2]*9)/sum_MZ
F_M2_Q2 = [F1_Vb2_Q2, F2_Vb2_Q2, F3_Vb2_Q2 ]

print('Base shear force for Mode1 = ' + str(F_M1_Q2))
print('Base shear force for Mode2 = ' + str(F_M2_Q2))

### for only one ritz vector
r1_r1 = [[1], [2], [3]]
r1_r1 = np.array(r1_r1)
## M_hat and K_hat 
M_hat_r1 = np.matmul(np.matmul(r1_r1.transpose(), M_matrix), r1_r1)
K_hat_r1 = np.matmul(np.matmul(r1_r1.transpose(), K_matrix), r1_r1)
## eig value problem
na_f_Q2_r1 = np.sqrt(K_hat_r1/M_hat_r1)
eig_V_Q2_r1 = r1_r1

period_Q2_r1 = 2*ma.pi/na_f_Q2_r1
T1_Q2_r1 = period_Q2_r1
V1_Q2_r1 = eig_V_Q2_r1

index_T1_Q2_r1 = Bsearch(left, right, T1_Q2_r1)


###  peak displacement for question 2
def cal_displacement_Q2_r1(dis, v, M):
	d = np.ones((3,1))
	m_bar = np.dot(np.dot(v.transpose(), M), v) # modal mass
	l = np.dot(np.dot(v.transpose(), M), d)/m_bar # load participation
	U = l*dis*v # displacement
	return U

peak_dis_1_Q2_r1 = find_peak_acc(index_T1_Q2_r1, dis_spectrum, T1_Q2_r1)

peak_dis_M1_Q2_r1 = cal_displacement_Q2(peak_dis_1_Q2_r1, V1_Q2_r1, M_matrix)
peak_dis_M1_Q2_r1 = np.array(peak_dis_M1_Q2_r1).reshape((1,3))
print('						')
print('Peak displacement vector of mode1(Q2--r1) = ' + str(peak_dis_M1_Q2_r1))

peak_acc_1_Q2_r1 = find_peak_acc(index_T1_Q2_r1, acc_spectrum, T1_Q2_r1)

sum_MZ = (M_matrix[0,0]*3 + M_matrix[1,1]*6 + M_matrix[2,2]*9)
sum_M = M_matrix[0,0]+M_matrix[1,1]+M_matrix[2,2]

V_b1_Q2_r1 = sum_M*peak_acc_1_Q2_r1

## base shear force for three modes
F1_Vb1_Q2_r1, F2_Vb1_Q2_r1, F3_Vb1_Q2_r1 = V_b1_Q2_r1 * (M_matrix[0,0]*3)/sum_MZ, V_b1_Q2_r1 * (M_matrix[1,1]*6)/sum_MZ, V_b1_Q2_r1 * (M_matrix[2,2]*9)/sum_MZ
F_M1_Q2_r1 = [F1_Vb1_Q2_r1, F2_Vb1_Q2_r1, F3_Vb1_Q2_r1]
F_M1_Q2_r1 = np.array(F_M1_Q2_r1).reshape((1,3))
print('Base shear force for Mode1(Q2--r1) = ' + str(F_M1_Q2_r1))
