# Project	: MathsAgainstTheCOVID19, Euler's method for epidemics
# Author    : Víctor Loureiro Sancho
# Created   : April, 2020

import numpy as np
import matplotlib.pyplot as plt
import csv

# Initializations
dt = 1/4  		# timestep Delta t
t_init = 0      # initial time
t_end = 35     	# stopping time

n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps

X = np.zeros(5)                   # create space for current X=[S,I,R,D]^T
dXdt = np.zeros(5)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((5,n_steps+1))   # create a storage array for X=[S,I,R,D]^T

S_init = np.linspace(205000,210000,5) 		# initial population of S
I_init = 6904           		  			# initial population of I
H_init = 3404           		  			# initial population of H
R_init = 517            		  			# initial population of R
D_init = 277            		  			# initial population of D

# SIHRD Model Parameters 
BETA = np.linspace(0.29,0.3,5) 	# contagion probability
GAMMA_I = np.linspace(0,0.1,5)		# (1/Duracion)
GAMMA_H = np.linspace(0.06,0.07,5)	# (1/Duracion)
ALFA_I = np.linspace(0.07,0.075,5)	# hospital probability
MU = np.linspace(0.02,0.025,5)		# lethality

# Calibration Process with known data
I_real = []; H_real = []; R_real = []; D_real = []

with open('../input/sihrd.csv') as File:  
	reader = csv.DictReader(File)
	for row in reader:
		fila = str(row).split(';')
		for i in range(len(fila)):
			fila[i] = fila[i].replace("'}","")
			fila[i] = fila[i].replace("{'","")
		I_real.append(int(fila[5]))
		H_real.append(int(fila[6]))
		R_real.append(int(fila[7]))
		D_real.append(int(fila[8]))

max_average_accuracy = 0; max_accuracy_i = 0; max_accuracy_h = 0; max_accuracy_r = 0; max_accuracy_d = 0
best_si = 0; best_beta = 0; best_gamma_i = 0; best_gamma_h = 0; best_alfa_i = 0; best_mu = 0

for Si in S_init:

	for b in BETA:
		
		for gi in GAMMA_I:

			for ai in ALFA_I:

				for gh in GAMMA_H:

					for m in MU:

						X_arr[0,0] = Si     		  # add the initial S to the storage array
						X_arr[1,0] = I_init     	  # add the initial I to the storage array
						X_arr[2,0] = H_init     	  # add the initial H to the storage array
						X_arr[3,0] = R_init     	  # add the initial R to the storage array
						X_arr[4,0] = D_init     	  # add the initial D to the storage array
						t_arr[0] = t_init       	  # add the initial t to the storage array

						# Euler's Method for SIR Model
						for i in range (1, n_steps + 1):
						    
						    t = t_arr[i-1]      # load the time
						    S = X_arr[0,i-1]    # load the value of S
						    I = X_arr[1,i-1]    # load the value of I
						    H = X_arr[2,i-1]    # load the value of H
						    R = X_arr[3,i-1]    # load the value of R
						    D = X_arr[4,i-1]    # load the value of D
						    N = S + I + H + R + D 

						    X[0] = S            # fill current state vector X=[S,I,H,R,D]^T
						    X[1] = I
						    X[2] = H
						    X[3] = R
						    X[4] = D

						    dSdt = -b*S*I/N   			# calculate the derivative dS/dt
						    dIdt = b*S*I/N - I*(gi+ai) 	# calculate the derivative dI/dt
						    dHdt = I*ai - H*(gh+m)  	# calculate the derivative dH/dt
						    dRdt = gh*H + gi*I 			# calculate the derivative dR/dt
						    dDdt = m*H     				# calculate the derivative dD/dt

						    dXdt[0] = dSdt      	# fill derivative vector dX/dt
						    dXdt[1] = dIdt
						    dXdt[2] = dHdt         
						    dXdt[3] = dRdt
						    dXdt[4] = dDdt

						    Xnew = X + dt*dXdt      # calculate X on next time step
						    X_arr[:,i] = Xnew       # store Xnew 
						    t_arr[i] = t + dt       # store new t-value 

						# Accuracy
						accuracy_i = 0; accuracy_h = 0; accuracy_r = 0; accuracy_d = 0;
						suma_i = 0;	suma_h = 0; suma_r = 0; suma_d = 0;
						
						for i in range(len(I_real)):
							suma_i += 1 - abs((I_real[i]-X_arr[1,int(i/dt)])/I_real[i])
							suma_h += 1 - abs((H_real[i]-X_arr[2,int(i/dt)])/H_real[i])
							suma_r += 1 - abs((R_real[i]-X_arr[3,int(i/dt)])/R_real[i])
							suma_d += 1 - abs((D_real[i]-X_arr[4,int(i/dt)])/D_real[i])

						accuracy_i = suma_i/len(I_real); accuracy_h = suma_h/len(H_real); accuracy_r = suma_r/len(R_real); accuracy_d = suma_d/len(D_real);
						average_accuracy = (accuracy_i+accuracy_h+accuracy_r+accuracy_d)/4

						if(average_accuracy>max_average_accuracy):
							
							max_accuracy_i = accuracy_i; max_accuracy_h = accuracy_h; max_accuracy_r = accuracy_r; max_accuracy_d = accuracy_d; max_average_accuracy = average_accuracy;
							best_si = Si; best_beta = b; best_gamma_i = gi; best_gamma_h = gh; best_alfa_i = ai; best_mu = m;
							print("[S(0) = " + str(best_si) + "\tBETA = " + str(best_beta) + "\tGAMMA_I = " + str(best_gamma_i) + "\tGAMMA_H = " + str(best_gamma_h) + "\tALFA_I = " + str(best_alfa_i) + "\tMU = " + str(best_mu) + "]\nI ACCURACY = " + str(max_accuracy_i) + "\nH ACCURACY = " + str(max_accuracy_h) + "\nR ACCURACY = " + str(max_accuracy_r) + "\nD ACCURACY = " + str(max_accuracy_d) + "\nAVERAGE ACCURACY = " + str(max_average_accuracy) + "\n")

# Euler's Calibrated Method for SIHRD Model
t_end = 150    	# stopping time
n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps

X = np.zeros(5)                   # create space for current X=[S,I,R,D]^T
dXdt = np.zeros(5)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((5,n_steps+1))   # create a storage array for X=[S,I,R,D]^T

S_init = best_si				  # initial population of S

# SIRD Model Parameters 
BETA = best_beta	 		# contagion probability
GAMMA_I = best_gamma_i		# (1/Duracion)
GAMMA_H = best_gamma_h		# (1/Duracion)
ALFA_I = best_alfa_i		# hospital probability
MU = best_mu				# lethality

X_arr[0,0] = S_init   		  # add the initial S to the storage array
X_arr[1,0] = I_init     	  # add the initial I to the storage array
X_arr[2,0] = H_init     	  # add the initial H to the storage array
X_arr[3,0] = R_init     	  # add the initial R to the storage array
X_arr[4,0] = D_init     	  # add the initial D to the storage array
t_arr[0] = t_init       	  # add the initial t to the storage array

for i in range (1, n_steps + 1):
    
    t = t_arr[i-1]      # load the time
    S = X_arr[0,i-1]    # load the value of S
    I = X_arr[1,i-1]    # load the value of I
    H = X_arr[2,i-1]    # load the value of H
    R = X_arr[3,i-1]    # load the value of R
    D = X_arr[4,i-1]    # load the value of D
    N = S + I + H + R + D

    X[0] = S            # fill current state vector X=[S,I,H,R,D]^T
    X[1] = I
    X[2] = H
    X[3] = R
    X[4] = D

    dSdt = -BETA*S*I/N 						# calculate the derivative dS/dt
    dIdt = BETA*S*I/N - I*(ALFA_I+GAMMA_I)  # calculate the derivative dI/dt
    dHdt = I*ALFA_I - H*(GAMMA_H+MU)  		# calculate the derivative dH/dt
    dRdt = GAMMA_I*I + GAMMA_H*H			# calculate the derivative dR/dt
    dDdt = MU*H     						# calculate the derivative dD/dt

    dXdt[0] = dSdt      	# fill derivative vector dX/dt
    dXdt[1] = dIdt
    dXdt[2] = dHdt         
    dXdt[3] = dRdt
    dXdt[4] = dDdt

    Xnew = X + dt*dXdt      # calculate X on next time step
    X_arr[:,i] = Xnew       # store Xnew 
    t_arr[i] = t + dt       # store new t-value

# Plot the results
fig = plt.figure()
plt.plot(t_arr, X_arr[0,:], linewidth = 4, label="S(t)")  	# plot S vs. time
plt.plot(t_arr, X_arr[1,:], linewidth = 4, label="I(t)")  	# plot I vs. time
plt.plot(t_arr, X_arr[2,:], linewidth = 4, label="H(t)")  	# plot H vs. time
plt.plot(t_arr, X_arr[3,:], linewidth = 4, label="R(t)")  	# plot R vs. time
plt.plot(t_arr, X_arr[4,:], linewidth = 4, label="D(t)")  	# plot D vs. time

plt.title('Eulers Method for SIHRD Model', fontsize = 20)   	# set title
plt.xlabel('t (days)', fontsize = 20)                   		# name of horizontal axis
plt.ylabel('S(t), I(t), H(t), R(t) and D(t)', fontsize = 20)    # name of vertical axis

plt.grid(which = 'both')                                				# set the grid
plt.xticks(fontsize = 15)                               				# adjust the fontsize
plt.yticks(fontsize = 15)                               				# adjust the fontsize
plt.axis([t_init, t_end, 0, 200000])                    				# set the range of the axes
plt.legend(fontsize=15)                                 				# show the legend
fig.savefig('../output/sihrd.png', dpi=fig.dpi, bbox_inches = "tight") 	# save the figure as .png
plt.show()                                              				# necessary for some platforms