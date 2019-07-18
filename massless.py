# We'll reproduce numerically the field equations (x-direction for the particle). Reference: Florian Scheck, Classical field theory on electrodynamics, Non - Abelian Gauge theories and gravitation (Section: Relativistic kinematics and dynamics, pg 123); F. Rohrlich, Classical charged particles (section: The Maxwell - Lorentz Field, pag 74)
import numpy as np
import random
import matplotlib.pyplot as plt
e = 1.0		#[c]
c = 1.0	#[m/s]
ti = 0.0		#Initial time of the observation
tf = 100.0		#Final time of the observation
N = 1000
Time = np.linspace(ti, tf, N)		#Respect to reference frame K which is at rest
#Function: Dot product
def dot(vector1, vector2):
	dot = -c**2.0*(vector1[0]*vector2[0])
	for i in range(len(vector1) - 1):
		dot += vector1[i+1]*vector2[i+1]
	return dot
#Function for the path of the particle. Insert here the function that you want to use for the path. We define the derivatives of the function
def path(k,delta,t):
	return k/(1.0 + (t/delta)**2.0)						#k*np.tanh(t/delta)
def derivative(k,delta,t):
	return -2.0*k/(delta**2.0)*t/((1.0 + (t/delta)**2.0)**2.0)				#k/(delta*np.cosh(t/delta)**2.0)
def second_derivative(k,delta,t):
	return -2.0*k/(delta**2.0)*(1.0 -3.0*(t/delta)**2.0)/((1.0 + (t/delta)**2.0)**3.0)	#-2.0*k/(delta)**2.0*np.tanh(t/delta)/(np.cosh(t/delta)**2.0)	
def gamma(deriva):			#Insert here the derivative of your path
	return 1.0/((1.0-(deriva/c)**2.0)**(0.5))
#Now, we establish the points where we want to calculate the electromagnetic field. Before, we establish the condition over the function's parameters: k/delta < c. We take k = c/2 and delta = 1
k = c/2.0
delta = 1.0
#We have the time coordinate. To establish the spatial coordinates we'll use a random function for N points inside a sphere of Radius R.
R = 20
field = [Time, np.zeros(N), np.zeros(N),np.zeros(N)]
for i in range(N):
	x1 = R	#random.uniform(0.0,R)
	x2 = R	#random.uniform(0.0,(R**2.0 - x1**2.0)**0.5)
	x3 = R	#random.uniform(0.0,(R**2.0 - x1**2.0 - x2**2.0)**0.5)
	field[1][i] = x1
	field[2][i] = x2
	field[3][i] = x3
#Now we'll obtain the zeros of the function R**2.0 = 0 (According to Scheck and Rohrlich)
particle = []
def function_zeros(field,t):
	return dot(field,field)-2.0*(-c**2.0*field[0]*t + field[1]*path(k,delta,t)) - c**2.0*t**2.0 + path(k,delta,t)**2.0
cota = -1000.0	#For the time Interval
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
null_vector = []
for i in range(4):
	null_vector.append(field[i] - particle[i])
#Proof R**2 = 0.0
print dot([null_vector[0][6],null_vector[1][6],null_vector[2][6],null_vector[3][6]], [null_vector[0][6],null_vector[1][6],null_vector[2][6],null_vector[3][6]])
#Now, we'll create the electromagnetic field tensor for the N points in the space-time. We'll use an array with N matrixes.
electromagnetic_field_N = []
for i in range(N):
	matrix = np.zeros((4,4))
	for j in range(4):
		for m in range(4):
			matrix[j][m] = e/(dot([particle_v[0][i],particle_v[1][i],particle_v[2][i],particle_v[3][i]], [null_vector[0][i],null_vector[1][i],null_vector[2][i],null_vector[3][i]])**3.0)*((null_vector[j][i]*particle_a[m][i] - null_vector[m][i]*particle_a[j][i])*dot([particle_v[0][i],particle_v[1][i],particle_v[2][i],particle_v[3][i]], [null_vector[0][i],null_vector[1][i],null_vector[2][i],null_vector[3][i]]) - (null_vector[j][i]*particle_v[m][i] - null_vector[m][i]*particle_v[j][i])*(c**2.0 + dot([particle_a[0][i],particle_a[1][i],particle_a[2][i],particle_a[3][i]], [null_vector[0][i],null_vector[1][i],null_vector[2][i],null_vector[3][i]])))
	electromagnetic_field_N.append(matrix)

Ex = np.zeros(N)
Ey = np.zeros(N)
Ez = np.zeros(N)
for i in range(N):
	Ex[i] = -electromagnetic_field_N[i][0,1]
	Ey[i] = -electromagnetic_field_N[i][0,2]
	Ez[i] = -electromagnetic_field_N[i][0,3]
Bx = np.zeros(N)
By = np.zeros(N)
Bz = np.zeros(N)
for i in range(N):
	Bx[i] = -electromagnetic_field_N[i][2,3]
	By[i] = electromagnetic_field_N[i][1,3]
	Bz[i] = -electromagnetic_field_N[i][1,2]
plt.figure()
plt.plot(Time,Ex,c ="r", label = "Ex")
plt.plot(Time,Ey,c ="b", label = "Ey")
plt.plot(Time,Ez,c ="g", label = "Ez")
plt.legend()
plt.savefig("Electric_Field_Acceleration.jpg")
plt.close()

plt.figure()
plt.plot(Time,Bx,c ="r", label = "Bx")
plt.plot(Time,By,c ="b", label = "By")
plt.plot(Time,Bz,c ="g", label = "Bz")
plt.legend()
plt.savefig("Magnetic_Field_Acceleration.jpg")
plt.close()


#Now, we'll verify that our code is right. We need the expressions for the electromagnetic field (Of corse, we use these in the other part in the definition of the field tensor)

normal_null_vector = np.zeros((N,3))
xn = np.zeros(N)
yn = np.zeros(N)
zn = np.zeros(N)
norm = np.zeros(N)
for i in range(N):
	norm[i] = ((field[1][i] - particle[1][i])**2.0 + (field[2][i] - particle[2][i])**2.0 + (field[3][i] - particle[3][i])**2.0)**(0.5)
	xn[i] = (field[1][i] - particle[1][i])/norm[i]
	yn[i] = (field[2][i] - particle[2][i])/norm[i]
	zn[i] = (field[3][i] - particle[3][i])/norm[i]
normal_null_vector[:,0] = xn
normal_null_vector[:,1] = yn
normal_null_vector[:,2] = zn

beta = np.zeros((N,3))			#Velocity
beta[:,0] = derivative(k,delta,zeros_time)/c

beta_dot = np.zeros((N,3))		#Acceleration
beta_dot[:,0] = second_derivative(k,delta,zeros_time)/c

gamma_value = gamma(derivative(k,delta,zeros_time))

electric_field_a = np.zeros((N,3))
for i in range(N):
	electric_field_a[i,:] = e/c*(np.cross(normal_null_vector[i,:],np.cross((normal_null_vector[i,:] - beta[i,:]), beta_dot[i,:])))/(norm[i]*(1-np.inner(beta[i,:],normal_null_vector[i,:]))**3.0)
#Now, we'll calculate the radiation emitted by an accelerated charge. For this part, we use the definitions above and the radial component of the pointying vector (See Classical Electrpdynamics, Jackson). We define the array radiation which contains the radiation emitted for the charged particle measured in a fix point in the space during a time T. 
radiation = np.zeros(N)
def norma(vector):
	value = 0.0
	for s in range(len(vector)):
		value += vector[s]**2.0
	return value**(0.5)
for i in range(N):
	radiation[i] = e**2.0/(4.0*np.pi*c)*1.0/(norm[i]**2.0)*1.0/(abs(1.0 - np.inner(beta[i,:],normal_null_vector[i,:]))**6.0)*norma(np.cross(normal_null_vector[i,:],np.cross((normal_null_vector[i,:] - beta[i,:]), beta_dot[i,:])))**2.0

plt.figure()
plt.plot(Time,radiation)
plt.show()
'''
electric_field = np.zeros((N,3))

for i in range(N):
	electric_field[i,:] = e*(normal_null_vector[i,:] - beta[i,:])/(gamma_value[i]**2.0*norm[i]**2.0*(1.0 - np.inner(beta[i,:],normal_null_vector[i,:]))) + e/c*(np.cross(normal_null_vector[i,:],np.cross((normal_null_vector[i,:] - beta[i,:]), beta_dot[i,:])))/(norm[i]*(1-np.inner(beta[i,:],normal_null_vector[i,:]))**3.0)

plt.figure()
plt.plot(Time,electric_field[:,0],c = "r", label = "Ex")
plt.plot(Time,electric_field[:,1],c = "b", label = "Ey")
plt.plot(Time,electric_field[:,2],c = "g", label = "Ez")
plt.legend()
plt.show()
plt.close()


magnetic_field = np.zeros((N,3))
for i in range(N):
	magnetic_field[i,:] = np.cross(normal_null_vector[i,:],electric_field[i,:])

plt.figure()
plt.plot(Time,magnetic_field[:,0],c = "r", label = "Bx")
plt.plot(Time,magnetic_field[:,1],c = "b", label = "By")
plt.plot(Time,magnetic_field[:,2],c = "g", label = "Bz")
plt.legend()
plt.show()
plt.close()
'''
	






