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

X = np.zeros(7)                   # create space for current X=[S,L,I,H,U,R,D]^T
dXdt = np.zeros(7)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((7,n_steps+1))   # create a storage array for X=[S,L,I,H,U,R,D]^T

# Known Data Acquisition
I_real = []; H_real = []; U_real = []; R_real = []; D_real = [];
with open('../input/slihurd.csv') as File:  
	reader = csv.DictReader(File)
	for row in reader:
		fila = str(row).split(';')
		for i in range(len(fila)):
			fila[i] = fila[i].replace("'}","")
			fila[i] = fila[i].replace("{'","")
		I_real.append(int(fila[6]))
		H_real.append(int(fila[7]))
		U_real.append(int(fila[8]))
		R_real.append(int(fila[9]))
		D_real.append(int(fila[10]))

S_init = np.linspace(275000,350000,10)	# initial population of S
L_init = np.linspace(0,50000,10)		# initial population of L
I_init = I_real[0]         		  		# initial population of I
H_init = H_real[0]         		  		# initial population of H
U_init = U_real[0]         		  		# initial population of U
R_init = R_real[0]         		  		# initial population of R
D_init = D_real[0]         		  		# initial population of D

# SLIHURD Model Parameters 
BETA = np.linspace(0.1,0.3,3)		# contagion probability
GAMMA_L = np.linspace(0.1,0.3,2)	# (1/Duracion)
GAMMA_I = np.linspace(0.005,0.01,2)	# (1/Duracion)
GAMMA_H = np.linspace(0.1,0.15,2)	# (1/Duracion)
GAMMA_U = np.linspace(0.01,0.05,2)	# (1/Duracion)
ALFA_I = np.linspace(0.05,0.1,2)	# Hospital probability
ALFA_H = np.linspace(0.02,0.07,2)	# ICUed probability
MU_H = np.linspace(0.005,0.01,2)	# hospitalized lethality
MU_U = np.linspace(0.05,0.1,2)		# ICUed lethality

# Calibration Process with known data
best_si = 0; best_li = 0; best_beta = 0; best_gamma_l = 0; best_gamma_i = 0; best_gamma_h = 0; best_gamma_u = 0; best_alfa_i = 0; best_alfa_h = 0; best_mu_h = 0; best_mu_u = 0;
max_accuracy_i = 0; max_accuracy_h = 0; max_accuracy_u = 0; max_accuracy_r = 0; max_accuracy_d = 0; max_average_accuracy = 0;

for Si in S_init:
	for Li in L_init:
		for b in BETA:
			for gl in GAMMA_L:
				for gi in GAMMA_I:
					for gh in GAMMA_H:
						for gu in GAMMA_U:
							for ai in ALFA_I:
								for ah in ALFA_H:
									for mh in MU_H:
										for mu in MU_U:

											X_arr[0,0] = Si     		# add the initial S to the storage array
											X_arr[1,0] = Li     	  	# add the initial L to the storage array
											X_arr[2,0] = I_init     	# add the initial I to the storage array
											X_arr[3,0] = H_init     	# add the initial H to the storage array
											X_arr[4,0] = U_init     	# add the initial U to the storage array
											X_arr[5,0] = R_init     	# add the initial R to the storage array
											X_arr[6,0] = D_init     	# add the initial D to the storage array
											t_arr[0] = t_init       	# add the initial t to the storage array

											# Euler's Method for SLIHURD Model
											for i in range (1, n_steps + 1):
											    
											    t = t_arr[i-1]      # load the time
											    S = X_arr[0,i-1]    # load the value of S
											    L = X_arr[1,i-1]    # load the value of L
											    I = X_arr[2,i-1]    # load the value of I
											    H = X_arr[3,i-1]    # load the value of H
											    U = X_arr[4,i-1]    # load the value of U
											    R = X_arr[5,i-1]    # load the value of R
											    D = X_arr[6,i-1]    # load the value of D
											    N = S + L + I + H + U + R + D

											    X[0] = S            # fill current state vector X=[S,L,I,H,U,R,D]^T
											    X[1] = L
											    X[2] = I
											    X[3] = H
											    X[4] = U
											    X[5] = R
											    X[6] = D

											    dSdt = -b*S*L/N   			# calculate the derivative dS/dt
											    dLdt = b*S*L/N - L*gl  		# calculate the derivative dL/dt
											    dIdt = L*gl - I*(gi+ai) 	# calculate the derivative dI/dt
											    dHdt = I*ai - H*(gh+ah+mh)	# calculate the derivative dH/dt
											    dUdt = H*ah - U*(gu+mu)		# calculate the derivative dU/dt
											    dRdt = I*gi + H*gh + U*gu 	# calculate the derivative dR/dt
											    dDdt = H*mh + U*mu 			# calculate the derivative dD/dt

											    dXdt[0] = dSdt      	# fill derivative vector dX/dt
											    dXdt[1] = dLdt
											    dXdt[2] = dIdt
											    dXdt[3] = dHdt
											    dXdt[4] = dUdt         
											    dXdt[5] = dRdt
											    dXdt[6] = dDdt

											    Xnew = X + dt*dXdt      # calculate X on next time step
											    X_arr[:,i] = Xnew       # store Xnew 
											    t_arr[i] = t + dt       # store new t-value 

											# Accuracy
											accuracy_i = 0; accuracy_h = 0; accuracy_u = 0; accuracy_r = 0; accuracy_d = 0;
											suma_i = 0;	suma_h = 0; suma_u = 0; suma_r = 0; suma_d = 0;
											
											for i in range(len(I_real)):
												suma_i += 1 - abs((I_real[i]-X_arr[2,int(i/dt)])/I_real[i])
												suma_h += 1 - abs((H_real[i]-X_arr[3,int(i/dt)])/H_real[i])
												suma_u += 1 - abs((U_real[i]-X_arr[4,int(i/dt)])/U_real[i])
												suma_r += 1 - abs((R_real[i]-X_arr[5,int(i/dt)])/R_real[i])
												suma_d += 1 - abs((D_real[i]-X_arr[6,int(i/dt)])/D_real[i])

											accuracy_i = suma_i/len(I_real); accuracy_h = suma_h/len(H_real); accuracy_u = suma_u/len(U_real); accuracy_r = suma_r/len(R_real); accuracy_d = suma_d/len(D_real);
											average_accuracy = (accuracy_i+accuracy_h+accuracy_u+accuracy_r+accuracy_d)/5

											if(average_accuracy>max_average_accuracy):
												
												max_accuracy_i = accuracy_i; max_accuracy_h = accuracy_h; max_accuracy_u = accuracy_u; max_accuracy_r = accuracy_r; max_accuracy_d = accuracy_d; max_average_accuracy = average_accuracy;
												best_si = Si; best_li = Li; best_beta = b; best_gamma_l = gl; best_gamma_i = gi; best_gamma_h = gh; best_gamma_u = gu; best_alfa_i = ai; best_alfa_h = ah; best_mu_h = mh; best_mu_u = mu;
												print("[S(0) = " + str(best_si) + "\tL(0) = " + str(best_li) + "\tBETA = " + str(best_beta) + "\tGAMMA_L = " + str(best_gamma_l) + "\tGAMMA_I = " + str(best_gamma_i) + "\tGAMMA_H = " + str(best_gamma_h) + "\tGAMMA_U = " + str(best_gamma_u) + "\tALFA_I = " + str(best_alfa_i) + "\tALFA_H = " + str(best_alfa_h) + "\tMU_H = " + str(best_mu_h) + "\tMU_U = " + str(best_mu_u) + "]\nI ACCURACY = " + str(max_accuracy_i) + "\nH ACCURACY = " + str(max_accuracy_h) + "\nU ACCURACY = " + str(max_accuracy_u) + "\nR ACCURACY = " + str(max_accuracy_r) + "\nD ACCURACY = " + str(max_accuracy_d) + "\nAVERAGE ACCURACY = " + str(max_average_accuracy) + "\n")

# Euler's Calibrated Method for SLIHURD Model
t_end = 150    							# stopping time
n_steps = int(round((t_end-t_init)/dt)) # total number of timesteps

X = np.zeros(7)                   # create space for current X=[S,L,I,H,U,R,D]^T
dXdt = np.zeros(7)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a storage array for t
X_arr = np.zeros((7,n_steps+1))   # create a storage array for X=[S,L,I,H,U,R,D]^T

# SLIRD Model Parameters 
BETA = best_beta	 		# contagion probability
GAMMA_L = best_gamma_l		# (1/Duracion)
GAMMA_I = best_gamma_i		# (1/Duracion)
GAMMA_H = best_gamma_h		# (1/Duracion)
GAMMA_U = best_gamma_u		# (1/Duracion)
ALFA_I = best_alfa_i		# Hospital probability
ALFA_H = best_alfa_h		# ICU probability
MU_H = best_mu_h			# lethality
MU_U = best_mu_u			# lethality

X_arr[0,0] = best_si   		  # add the initial S to the storage array
X_arr[1,0] = best_li     	  # add the initial L to the storage array
X_arr[2,0] = I_init     	  # add the initial I to the storage array
X_arr[3,0] = H_init     	  # add the initial H to the storage array
X_arr[4,0] = U_init     	  # add the initial U to the storage array
X_arr[5,0] = R_init     	  # add the initial R to the storage array
X_arr[6,0] = D_init     	  # add the initial D to the storage array
t_arr[0] = t_init       	  # add the initial t to the storage array

for i in range (1, n_steps + 1):
    
    t = t_arr[i-1]      # load the time
    S = X_arr[0,i-1]    # load the value of S
    L = X_arr[1,i-1]    # load the value of L
    I = X_arr[2,i-1]    # load the value of I
    H = X_arr[3,i-1]    # load the value of H
    U = X_arr[4,i-1]    # load the value of U
    R = X_arr[5,i-1]    # load the value of R
    D = X_arr[6,i-1]    # load the value of D
    N = S + L + I + H + U + R + D

    X[0] = S            # fill current state vector X=[S,L,I,H,U,R,D]^T
    X[1] = L
    X[2] = I
    X[3] = H
    X[4] = U
    X[5] = R
    X[6] = D

    dSdt = -BETA*S*L/N   						# calculate the derivative dS/dt
    dLdt = BETA*S*L/N - L*GAMMA_L  				# calculate the derivative dL/dt
    dIdt = L*GAMMA_L - I*(GAMMA_I+ALFA_I) 		# calculate the derivative dI/dt
    dHdt = I*ALFA_I - H*(GAMMA_H+ALFA_H+MU_H)	# calculate the derivative dH/dt
    dUdt = H*ALFA_H - U*(GAMMA_U+MU_U)			# calculate the derivative dU/dt
    dRdt = I*GAMMA_I + H*GAMMA_H + U*GAMMA_U 	# calculate the derivative dR/dt
    dDdt = H*MU_H + U*MU_U 						# calculate the derivative dD/dt

    dXdt[0] = dSdt      	# fill derivative vector dX/dt
    dXdt[1] = dLdt
    dXdt[2] = dIdt         
    dXdt[3] = dHdt
    dXdt[4] = dUdt
    dXdt[5] = dRdt
    dXdt[6] = dDdt

    Xnew = X + dt*dXdt      # calculate X on next time step
    X_arr[:,i] = Xnew       # store Xnew 
    t_arr[i] = t + dt       # store new t-value

# Plot the results
fig = plt.figure()
plt.plot(t_arr, X_arr[0,:], linewidth = 4, label="S(t)")  	# plot S vs. time
plt.plot(t_arr, X_arr[1,:], linewidth = 4, label="L(t)")  	# plot L vs. time
plt.plot(t_arr, X_arr[2,:], linewidth = 4, label="I(t)")  	# plot I vs. time
plt.plot(t_arr, X_arr[3,:], linewidth = 4, label="H(t)")  	# plot H vs. time
plt.plot(t_arr, X_arr[4,:], linewidth = 4, label="U(t)")  	# plot U vs. time
plt.plot(t_arr, X_arr[5,:], linewidth = 4, label="R(t)")  	# plot R vs. time
plt.plot(t_arr, X_arr[6,:], linewidth = 4, label="D(t)")  	# plot D vs. time

plt.title('Eulers Method for SLIHURD Model', fontsize = 20)   				# set title
plt.xlabel('t (days)', fontsize = 20)                   					# name of horizontal axis
plt.ylabel('S(t), L(t), I(t), H(t), U(t), R(t) and D(t)', fontsize = 20)    # name of vertical axis

plt.grid(which = 'both')                                					# set the grid
plt.xticks(fontsize = 15)                               					# adjust the fontsize
plt.yticks(fontsize = 15)                               					# adjust the fontsize
plt.axis([t_init, t_end, 0, 300000])                    					# set the range of the axes
plt.legend(fontsize=15)                                 					# show the legend
fig.savefig('../output/slihurd.png', dpi=fig.dpi, bbox_inches = "tight") 	# save the figure as .png
plt.show()                                              					# necessary for some platforms