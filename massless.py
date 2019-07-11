# We'll reproduce numerically the field equations (x-direction for the particle). Reference: Florian Scheck, Classical field theory on electrodynamics, Non - Abelian Gauge theories and gravitation (Section: Relativistic kinematics and dynamics, pg 123); F. Rohrlich, Classical charged particles (section: The Maxwell - Lorentz Field, pag 74)
import numpy as np
import random
import matplotlib.pyplot as plt
e = 1.0			#[c]
c = 1.0			#[m/s]
ti = 0.0		#Initial time of the observation
tf = 10.0		#Final time of the observation
N = 20
Time = np.linspace(ti, tf, N)		#Respect to reference frame K which is at rest
#Function: Dot product
def dot(vector1, vector2):
	dot = -1.0*(vector1[0]*vector2[0])
	i = 1
	for i in range(len(vector1) - 1):
		dot += vector1[i+1]*vector2[i+1]
	return dot
#Function for the path of the particle. Insert here the function that you want to use for the path. We define the derivatives of the function
def path(k,delta,t):
	return k*np.tanh(t/delta)
def derivative(k,delta,t):
	return k/(delta*np.cosh(t/delta))
def second_derivative(k,delta,t):
	return -k/(delta)**2.0*np.tanh(t/delta)/(np.cosh(t/delta))	
def gamma(deriva):			#Insert here the derivative of your path
	return 1.0/((1.0-(deriva/c)**2.0)**(0.5))

#Now, we establish the points where we want to calculate the electromagnetic field. Before, we establish the condition over the function's parameters: k/delta < c. We take k = c/2 and delta = 1
k = c/2.0
delta = 1.0
#We have the time coordinate. To establish the spatial coordinates we'll use a random function for N points inside a sphere of Radius R.
R = 10.0
field = [Time, np.zeros(N), np.zeros(N),np.zeros(N)]
for i in range(N):
	x1 = random.uniform(0.0,R)
	x2 = random.uniform(0.0,(R**2.0 - x1**2.0)**0.5)
	x3 = random.uniform(0.0,(R**2.0 - x1**2.0 - x2**2.0)**0.5)
	field[1][i] = x1
	field[2][i] = x2
	field[3][i] = x3
print field
#Now we define the orthogonal space-like vector v
orthogonal = [np.zeros(N),np.zeros(N),np.ones(N),np.zeros(N)]
#Now we'll obtain the zeros of the function R**2.0 = 0 (According to Scheck and Rohrlich)
particle = []
def function_zeros(field,t):
	return dot(field,field)-2.0*(-c**2.0*field[0]*t + field[1]*path(k,delta,t)) - c**2.0*t**2.0 + path(k,delta,t)**2.0
cota = -100.0	#For the time Interval
zeros_time = np.zeros(N)
for i in range(N):			#Zeros: Array with the information: Time, position particle
	time = random.uniform(cota + field[0][i], field[0][i])
	value_o = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],cota + field[0][i])
	value_i = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],time)
	if (abs(value_o) <= 0.01):
		zeros_time[i] = cota + field[0][i]	
	else:
		while(((value_i > 0.0 and value_o > 0.01) or (value_i < 0.0 and value_o < -0.01))):
			time = random.uniform(cota + field[0][i], field[0][i])
			value_i = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],time)
		if(abs(value_i) <= 0.01):
			zeros_time[i] = time
		else:			#Biyection method
			medium = (cota + field[0][i] + time)/2.0
			value_m = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],medium)
			while (abs(value_m) > 0.01):
				if((value_m < -0.01 and value_i > 0.01) or ((value_m > 0.01 and value_i < -0.01))):
					value_o = value_m
					medium = (medium + time)/2.0
					value_m = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],medium)
				elif((value_m < -0.01 and value_o > 0.01) or ((value_m > 0.01 and value_o < -0.01))):
					value_i = value_m
					time = medium
					medium = (medium + cota + field[0][i])/2.0
					value_m = function_zeros([field[0][i],field[1][i],field[2][i],field[3][i]],medium)
			zeros_time[i] = medium
#We obtaine the position for the particle such that field - particle = R is a null vector.
particle.append(zeros_time)		#Times where R is a null vector
positions = path(k,delta,zeros_time)
particle.append(positions)
particle.append(np.zeros(N))
particle.append(np.zeros(N))

#Finally, we want to calculate the field using that points. First, we need to calculate the acceleration, the velocity and rho
particle_v = gamma(derivative(k,delta,particle[0]))*[c*np.ones(N),derivative(k,delta,particle[0]),np.zeros(N),np.zeros(N)]
particle_a = gamma(derivative(k,delta,particle[0]))**2.0*[gamma(derivative(k,delta,particle[0]))**2.0/c*derivative(k,delta,particle[0])*second_derivative(k,delta,particle[0]),gamma(derivative(k,delta,particle[0]))**2.0/(c**2.0)*derivative(k,delta,particle[0])**2.0*second_derivative(k,delta,particle[0]) + second_derivative(k,delta,particle[0]),np.zeros(N),np.zeros(N)]

#Now, we calculate the null vector R
field[0] = c*field[0]
particle[0] =c*particle[0]
null_vector = []
for i in range(4):
	null_vector.append(field[i] - particle[i])
'''#Proof R**2 = 0.0
print dot([null_vector[0][6],null_vector[1][6],null_vector[2][6],null_vector[3][6]], [null_vector[0][6],null_vector[1][6],null_vector[2][6],null_vector[3][6]])'''
rho = np.zeros(N)
for i in range(N):
	rho[i] = dot([orthogonal[0][i],orthogonal[1][i],orthogonal[2][i],orthogonal[3][i]],[null_vector[0][i],null_vector[1][i],null_vector[2][i],null_vector[3][i]])
scalar_acceleration = np.zeros(N)
for i in range(N):
	scalar_acceleration[i] = dot([orthogonal[0][i],orthogonal[1][i],orthogonal[2][i],orthogonal[3][i]],[particle_a[0][i],particle_a[1][i],particle_a[2][i],particle_a[3][i]])
#Now, we'll create the electromagnetic field tensor for the N points in the space-time. We'll use an array with N matrixes.
electromagnetic_field_N = []
for i in range(N):
	matrix = np.zeros((4,4))
	for j in range(4):
		for k in range(4):
			matrix[j][k] = e/(c*(rho[i])**2.0)*(particle_v[j][i]*orthogonal[k][i] - particle_v[k][i]*orthogonal[j][i]) + e/(c**2.0*rho[i])*((particle_a[j][i]*particle_v[k][i] - particle_a[k][i]*particle_v[j][i])/c - orthogonal[j][i]*(scalar_acceleration[i]*particle_v[k][i]/c + particle_a[k][i]) + orthogonal[k][i]*(scalar_acceleration[i]*particle_v[j][i]/c + particle_a[j][i]))
	electromagnetic_field_N.append(matrix)














	

