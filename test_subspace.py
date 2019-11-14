"""
This is the test code for the subspace identification method.

"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

from functionsSID import estimateMarkovParameters
from functionsSID import estimateModel
from functionsSID import systemSimulate
from functionsSID import estimateInitial
from functionsSID import modelError

###############################################################################
# Define the model

# masses, spring and damper constants
m1=20  ; m2=20   ; k1=1000  ; k2=2000 ; d1=1  ; d2=5; 
# define the continuous-time system matrices
Ac=np.matrix([[0, 1, 0, 0],[-(k1+k2)/m1 ,  -(d1+d2)/m1 , k2/m1 , d2/m1 ], [0 , 0 ,  0 , 1], [k2/m2,  d2/m2, -k2/m2, -d2/m2]])
Bc=np.matrix([[0],[0],[0],[1/m2]])
Cc=np.matrix([[1, 0, 0, 0]])
# end of model definition
###############################################################################

###############################################################################
# parameter definition

r=1; m=1 # number of inputs and outputs
# total number of time samples
time=300
# discretization constant
sampling=0.05

# model discretization
I=np.identity(Ac.shape[0]) # this is an identity matrix
A=inv(I-sampling*Ac)
B=A*sampling*Bc
C=Cc

# check the eigenvalues
eigen_A=np.linalg.eig(Ac)[0]
eigen_Aid=np.linalg.eig(A)[0]

# define an input sequence and initial state for the identification
input_ident=np.random.rand(1,time)
x0_ident=np.random.rand(4,1)

#define an input sequence and initial state for the validation
input_val=np.random.rand(1,time)
x0_val=np.random.rand(4,1)

# simulate the discrete-time system to obtain the input-output data for identification and validation
Y_ident, X_ident=systemSimulate(A,B,C,input_ident,x0_ident)
Y_val, X_val=systemSimulate(A,B,C,input_val,x0_val)

#  end of parameter definition
###############################################################################

###############################################################################
# model estimation and validation

# estimate the Markov parameters
past_value=10 # this is the past window - p 
Markov,Z, Y_p_p_l =estimateMarkovParameters(input_ident,Y_ident,past_value)

# estimate the system matrices
model_order=3 # this is the model order \hat{n}
Aid,Atilde,Bid,Kid,Cid,s_singular,X_p_p_l = estimateModel(input_ident,Y_ident,Markov,Z,past_value,past_value,model_order)  

plt.plot(s_singular, 'x',markersize=8)
plt.xlabel('Singular value index')
plt.ylabel('Singular value magnitude')
plt.yscale('log')
#plt.savefig('singular_values.png')
plt.show()

# estimate the initial state of the validation data
h=10 # window for estimating the initial state
x0est=estimateInitial(Aid,Bid,Cid,input_val,Y_val,h)

# simulate the open loop model 
Y_val_prediction,X_val_prediction = systemSimulate(Aid,Bid,Cid,input_val,x0est)

# compute the errors
relative_error_percentage, vaf_error_percentage, Akaike_error = modelError(Y_val,Y_val_prediction,r,m,30)
print('Final model relative error %f and VAF value %f' %(relative_error_percentage, vaf_error_percentage))

# plot the prediction and the real output 
plt.plot(Y_val[0,:100],'k',label='Real output')
plt.plot(Y_val_prediction[0,:100],'r',label='Prediction')
plt.legend()
plt.xlabel('Time steps')
plt.ylabel('Predicted and real outputs')
#plt.savefig('results3.png')
plt.show()

#               end of code
###############################################################################