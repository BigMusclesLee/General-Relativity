import numpy as np
import matplotlib.pyplot as plt
plt.clf() # Clears figure

# parameters
R0 = 12 # Initial radius, has to be greater than or equal to 6
alpha = 16 # Angular momentum, has to be greater than 12, for circular orbit set to R0/(1-3/R0)
dRdTau0 = 0 # Boundary condition, no initial radial motion at tau=0
M = 100 # Mass

def d2RdTau2(R): # Defines the derivative as a function
    return -1/((R**2)*(M**2)) + alpha/((R**3)*(M**2)) - 3*alpha/((R**4)*(M**2))

endTau = 5000 # End proper time
steps = 20000 # Steps for controllable degree of accuracy
taus = np.linspace(0, endTau, steps) # An array with size of steps with equally spaces numbers between 0 and the end tau
h = taus[1]-taus[0] # Step length

R = np.zeros(len(taus)) # Creates an array of 0's with size of how many taus there are
dRdTau = np.zeros(len(taus)) # Same as above but for radial motion
R[0] = R0 # Sets initial radius
dRdTau[0] = dRdTau0 # Sets initila radial motion

for i in range(len(taus)-1): # Start of Runge-Kutta 4
    k1_R = dRdTau[i] # Setting initial radial motion
    k1_dRdTau = d2RdTau2(R[i]) # Setting initial radial acceleration
    k2_R = dRdTau[i] + 0.5*h*k1_dRdTau # Estimated radial motion halfway to next step
    k2_dRdTau = d2RdTau2(R[i] + 0.5*h*k1_R) # Estimated radial acceleration halfway to next step
    k3_R = dRdTau[i] + 0.5*h*k2_dRdTau # Use previous estimate to make another
    k3_dRdTau = d2RdTau2(R[i] + 0.5*0.5*h*k2_R) # Same as above
    k4_R = dRdTau[i] + h*k3_dRdTau # Estimated radial motion for end of step
    k4_dRdTau = d2RdTau2(R[i] + h*k3_R) # Estimated radial acceleration for end of step
    R[i+1] = R[i] + (h/6.0)*(k1_R + 2*k2_R + 2*k3_R +k4_R) # By combining all estimates we get a weighted average for R after the step
    dRdTau[i+1] = dRdTau[i] + (h/6.0)*(k1_dRdTau + 2*k2_dRdTau + 2*k3_dRdTau + k4_dRdTau) # By combining all esimates we get a weighted average for radial motion after the step

# Next section is for part 6 only:
x = np.zeros(len(taus)) 
y = np.zeros(len(taus))
for i in range(len(taus)):
    x[i] = R[i]*np.cos(0.1*taus[i]) # Uses given theta to get x
    y[i] = R[i]*np.sin(0.1*taus[i]) # Uses given theta to get y

plt.plot(x, y, label=fr"$x(\tau),\, y(\tau), \alpha={alpha},\, R_0={R0},\, M={M}$") # Plots parameterised x by y and labels parameters
plt.xlabel(r"x") # Labels x axis
plt.ylabel(r"y") # Labels y axis

plt.axhline(y=0, color='k') # Making axis clearer
plt.axvline(x=0,color='k') # Making axis clearer

plt.title("Radial Motion of a Particle Orbiting a Black Hole in Schwarzschild Spacetime") # Descriptive title
plt.grid(True) # Adds grid
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.15), fontsize=14) # Moving legend outside of graph and making it bigger
plt.tight_layout() # Makes graph fit nicer
figNum="Figure 16" # Figure number to label graph
plt.figtext(0.5,0.01,figNum,wrap=True,horizontalalignment='center',fontsize=14) # Adds the label with figure number
plt.savefig("plot.png") # Saves figure