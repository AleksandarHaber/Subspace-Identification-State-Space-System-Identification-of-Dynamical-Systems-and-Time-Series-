## Subspace Identification of Dynamical Systems and Time Series

**IMPORTANT NOTE: First, thoroughly read the license in the file called LICENSE.md!**

A detailed post accompanying the codes can be found on this webpage:
https://aleksandarhaber.com/introduction-to-subspace-system-identification-system-identification-tutorial/

Briefly speaking, the coded subspace identification method estimates a state-space model

x_{k+1}=Ax_{k}+Bu_{k}


y_{k}=Cx_{k}

or the Kalman innovation state space model

x_{k+1}=Ax_{k}+Bu_{k}+Ke_{k}


y_{k}=Cx_{k}+e_{k}

using only the input-output data sequences (u_{k},y_{k}) for k=0,1,2,..., N. Also, the method can estimate the model order. 

I have implemented a version of the Predictor Based Subspace IDentification (PBSID) method. The description of the algorithm can be found in the following papers:

A. Haber, Subspace Identification of Temperature Dynamics, https://arxiv.org/abs/1908.02379
https://arxiv.org/abs/1908.02379

If you want to refer to this implementation in your paper or report, cite that paper and this GitHub repository.
Explanation of the files:

"test_subspace.py" - is a driver code for testing the subspace identification method, you should start from here. 

"discretization_test.py" - is used to test the backward Euler method for discretizing a dynamical system

"functionsSID.py" - contains the functions used to estimate the model, the functions from this file are imported in "test_subspace.py" 

All the questions and comments should be addressed to me:

Aleksandar Haber 
aleksandar.haber@gmail.com

Also, if you are interested in improving the codes, please let me know. This implementation is not optimized for large-scale problems. 

