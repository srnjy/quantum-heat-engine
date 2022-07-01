''' Objective: Estimate the efficiency (eta) of 
     a 2-site tight-binding Hamiltonian '''
##################################################
# Note: The code works for interacting Hubbard Hamiltonian
# To make it non-interacting tight-binding, set interaction U = 0. 

from matplotlib.widgets import Slider  # import the Slider widget
import matplotlib.pyplot as plt
from math import pi
import numpy as np
from numpy import linalg as lin
import math

eps = 0.5 # orbital energy
U = 0  # interaction 
T_high = 3.0 
T_low = 1.0
dim = 16
H = np.zeros([dim, dim], float)
N_state = 16
flag_eigen = 0
fix = t1 = 1.0  # fixed parameter
var = t2 = 2.0  # varying parameter
print('var,fix=',var,fix)

# Initialize
E = np.zeros(N_state)
P = np.zeros(N_state)
P_high = np.zeros(N_state)
P_low = np.zeros(N_state)

# Function: Eigenenergy by numerically solving Hamiltonian matrix 
def Eigen_num(t):

  # Construct Hamiltonian matrix 
  H[1][1] = H[2][2] = H[3][3] = H[4][4] = eps
  H[1][2] = H[2][1] = H[3][4] = H[4][3] = -t

  H[5][5] = H[7][7] = H[8][8] = H[10][10] = H[10][10] = 2 * eps
  H[6][6] = H[9][9] = 2 * eps + U
  H[6][7] = H[7][6] = H[7][9] = H[9][7] = -t
  H[6][8] = H[8][6] = H[8][9] = H[9][8] = t

  H[11][11] = H[12][12] = H[13][13] = H[14][14] = 3 * eps + U
  H[11][12] = H[12][11] = H[13][14] = H[14][13] = t

  H[15][15] = 4 * eps + U

  E, v1 = lin.eig(H)

  return E # should be a list of 16 elements

# Function: Eigenenergy found analytically
def Eigen_ana(t):
   E[0] = E[5] = E[6] = E[10] = 0
   E[1] = E[3] = -t
   E[2] = E[4] = t
   E[7] = U
   E[8] = 0.5*(U-np.sqrt(U**2+16*t**2))
   E[9] = 0.5*(U+np.sqrt(U**2+16*t**2))
   E[11] = E[13] = U - t
   E[12] = E[14] = U + t
   E[15] = 2 * U

   return E

# Function: Partition Function
def PF(temp,variable):
        beta=1.0/temp
        print('CALLING EIGEN FROM PF, flag_eigen =',flag_eigen)
        if flag_eigen == 0:
            E = Eigen_ana(variable)
        else:
            E = Eigen_num(variable)

        Z = 0; N_lim = N_state
        for i in range(N_lim):
          Z += np.exp(-beta*E[i])
        print('Partition fn.=',Z)
        return Z

# Function: Probabilty array
def P_array(temp,variable): # fixvar => fixed variable

        beta=1.0/temp
        P = np.zeros(N_state)  # This initialization is essential
        Z = PF(temp,variable) # This finds E-values as well,
                                           #     check PF function.
        for j in range(N_state):
             P[j] = np.exp(-beta*E[j])/Z
             print('j,P[j],E[j] =', j,P[j],E[j]) 
        return P

# Function: Heat energy
def Heat(variable):
     print('variable=',variable)
     if flag_eigen == 0:
            E = Eigen_ana(variable)
     else:
            E = Eigen_num(variable)
     print('T_high,var,T_low,fix=',T_high,var,T_low,fix)
     Q = 0
     for i in range(N_state):
            Q += E[i]*(P_high[i] - P_low[i])
            print('Q=',Q)
     return Q


# ----------------------------------------- #
#         Main part: 
# ----------------------------------------- #
var_max = 6.0  # maximum value of variable
dvar = 0.1 # variable increment
Npoint = int((var_max-var)/dvar) # no of grid points

print( 'Finding probability for fixed t2 ...' )
print('T_low,fix=',T_low,fix)
P_low = P_array(T_low,fix)


f0 = open('eta_vs_var.dat', 'w')
for iter in range(Npoint):   # loop for the varying parameter var
          P_high = P_array(T_high,var) # updates in every itern
          Q_high = Heat(var)
          Q_low =  Heat(fix) 
          Work = Q_high - Q_low
          print('Q_high,Q_low,Work=',Q_high,Q_low,Work)
          #if (Work<0.0): 
           #   Work=0.0 
          effy = Work/Q_high
          print('var,effy =',var,effy)
          print ('{:4.3f} \t {:4.3f}'.format(var,effy), file=f0)
          #plt.plot(var,effy,c='b', linestyle='dashed')
          if iter == 0:
            plt.scatter(var,effy, c='r', marker='+', label='numerical')
          else:
            plt.scatter(var,effy, c='r', marker='+')
       
          #effy_ana = 1 - (fix/var)**2  # anayltical efficiency
          effy_ana = 1 - (fix/var)  # anayltical efficiency
          if iter == 0:
            plt.scatter(var,effy_ana, c='b', marker='x', label='$1-t_2/t_1$')
          else:
            plt.scatter(var,effy_ana, c='b', marker='x')
      
          var += dvar  # update variable
f0.close()

# Annotate
xannote=3.0; yannote=.2
bbox_props = dict(boxstyle="round", fc="pink", ec="0.5", alpha=0.9)
plt.text(xannote, yannote, "$t_2$ = {}, $U$ = {}, $T_H$ = {}, $T_L$ = {}".format(t2,U,T_high,T_low), size=8, bbox=bbox_props)


plt.title('Efficiency of a 2-site heat engine')
plt.grid(True)
plt.legend(loc='best')
plt.xlabel('$t_1$', size=18)
plt.ylabel('$\eta$', size=18)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.16,  top=0.93) # Adjust layout
plt.savefig("eta_vs_t1_2site.png")
plt.show()
