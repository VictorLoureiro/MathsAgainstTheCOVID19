# Project	: MathsAgainstTheCOVID19, Euler's method for epidemics
# Author    : Víctor Loureiro Sancho
# Created   : April, 2020

import numpy as np
import matplotlib.pyplot as plt
import csv

# Initializations
dt = 1/4   		# timestep Delta t
t_init = 0      # initial time
t_end = 35     	# stopping time

n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps

X = np.zeros(3)                   # create space for current X=[S,I,R]^T
dXdt = np.zeros(3)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((3,n_steps+1))   # create a storage array for X=[S,I,R]^T

S_init = np.linspace(160000,161000,10)  	# initial population of S
I_init = 6904           		  			# initial population of I
R_init = 517            		  			# initial population of R

# SIR Model Parameters 
BETA = np.linspace(0.23,0.24,10)	 	# contagion probability
GAMMA = np.linspace(0.026,0.027,10)	# (1/Duracion)

# Calibration Process with known data
I_real = [] 
R_real = []

with open('../input/sir.csv') as File:  
	reader = csv.DictReader(File)
	for row in reader:
		fila = str(row).split(';')
		I_real.append(int(fila[4]))
		R_real.append(int(fila[5]))

max_accuracy_i = 0; max_accuracy_r = 0; max_average_accuracy = 0
best_si = 0; best_beta = 0; best_gamma = 0

for Si in S_init:

	for b in BETA:
		
		for g in GAMMA:

			X_arr[0,0] = Si     		  # add the initial S to the storage array
			X_arr[1,0] = I_init     	  # add the initial I to the storage array
			X_arr[2,0] = R_init     	  # add the initial R to the storage array
			t_arr[0] = t_init       	  # add the initial t to the storage array

			# Euler's Method for SIR Model
			for i in range (1, n_steps + 1):
			    
			    t = t_arr[i-1]      # load the time
			    S = X_arr[0,i-1]    # load the value of S
			    I = X_arr[1,i-1]    # load the value of I
			    R = X_arr[2,i-1]    # load the value of R
			    N = S + I + R

			    X[0] = S            # fill current state vector X=[S,I,R]^T
			    X[1] = I
			    X[2] = R

			    dSdt = -b*S*I/N   		# calculate the derivative dS/dt
			    dIdt = b*S*I/N - I*g  	# calculate the derivative dI/dt
			    dRdt = g*I     			# calculate the derivative dR/dt

			    dXdt[0] = dSdt      	# fill derivative vector dX/dt
			    dXdt[1] = dIdt         
			    dXdt[2] = dRdt

			    Xnew = X + dt*dXdt      # calculate X on next time step
			    X_arr[:,i] = Xnew       # store Xnew 
			    t_arr[i] = t + dt       # store new t-value 

			# Accuracy
			accuracy_i = 0; accuracy_r = 0;
			suma_i = 0;	suma_r = 0;
			
			for i in range(len(I_real)):
				suma_i += 1 - abs((I_real[i]-X_arr[1,int(i/dt)])/I_real[i])

			for i in range(len(R_real)):
				suma_r += 1 - abs((R_real[i]-X_arr[2,int(i/dt)])/R_real[i])

			accuracy_i = suma_i/len(I_real)
			accuracy_r = suma_r/len(R_real)
			average_accuracy = (accuracy_i+accuracy_r)/2

			if(average_accuracy>max_average_accuracy):
				
				max_accuracy_i = accuracy_i
				max_accuracy_r = accuracy_r
				max_average_accuracy = average_accuracy
				best_si = Si
				best_beta = b
				best_gamma = g
				print("[S(0) = " + str(best_si) + "\tBETA = " + str(best_beta) + "\tGAMMA = " + str(best_gamma) + "]\nI ACCURACY = " + str(max_accuracy_i) + "\nR ACCURACY = " + str(max_accuracy_r) + "\nAVERAGE ACCURACY = " + str(max_average_accuracy) + "\n")

# Euler's Calibrated Method for SIR Model
t_end = 150    	# stopping time
n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps

X = np.zeros(3)                   # create space for current X=[S,I,R]^T
dXdt = np.zeros(3)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((3,n_steps+1))   # create a storage array for X=[S,I,R]^T

S_init = best_si				  # initial population of S

# SIR Model Parameters 
BETA = best_beta	 			# contagion probability
GAMMA = best_gamma				# (1/Duracion)

X_arr[0,0] = S_init   		  # add the initial S to the storage array
X_arr[1,0] = I_init     	  # add the initial I to the storage array
X_arr[2,0] = R_init     	  # add the initial R to the storage array
t_arr[0] = t_init       	  # add the initial t to the storage array

for i in range (1, n_steps + 1):
    
    t = t_arr[i-1]      # load the time
    S = X_arr[0,i-1]    # load the value of S
    I = X_arr[1,i-1]    # load the value of I
    R = X_arr[2,i-1]    # load the value of R
    N = S + I + R

    X[0] = S            # fill current state vector X=[S,I,R]^T
    X[1] = I
    X[2] = R

    dSdt = -BETA*S*I/N   			# calculate the derivative dS/dt
    dIdt = BETA*S*I/N - I*GAMMA  	# calculate the derivative dI/dt
    dRdt = GAMMA*I     				# calculate the derivative dR/dt

    dXdt[0] = dSdt      	# fill derivative vector dX/dt
    dXdt[1] = dIdt         
    dXdt[2] = dRdt

    Xnew = X + dt*dXdt      # calculate X on next time step
    X_arr[:,i] = Xnew       # store Xnew 
    t_arr[i] = t + dt       # store new t-value

# Plot the results
fig = plt.figure()
plt.plot(t_arr, X_arr[0,:], linewidth = 4, label="S(t)")  	# plot S vs. time
plt.plot(t_arr, X_arr[1,:], linewidth = 4, label="I(t)")  	# plot I vs. time
plt.plot(t_arr, X_arr[2,:], linewidth = 4, label="R(t)")  	# plot R vs. time

plt.title('Eulers Method for SIR Model', fontsize = 20)   	# set title
plt.xlabel('t (days)', fontsize = 20)                   	# name of horizontal axis
plt.ylabel('S(t), I(t) and R(t)', fontsize = 20)        	# name of vertical axis

plt.grid(which = 'both')                                				# set the grid
plt.xticks(fontsize = 15)                               				# adjust the fontsize
plt.yticks(fontsize = 15)                               				# adjust the fontsize
plt.axis([t_init, t_end, 0, 160000])                    				# set the range of the axes
plt.legend(fontsize=15)                                 				# show the legend
fig.savefig('../output/sir.png', dpi=fig.dpi, bbox_inches = "tight") 	# save the figure as .png
plt.show()                                              				# necessary for some platforms