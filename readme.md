## Subspace Identification of Dynamical Systems and Time Series

A detailed post accompanying the codes can be found on my GutHub webpage:
https://aleksandarhaber.github.io/machine_learning/2019/11/13/subspace-identification.html

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

Houtzager, I., van Wingerden, J. W., & Verhaegen, M. (2009, December). VARMAX-based closed-loop subspace model identification. 
In Proceedings of the 48h IEEE Conference on Decision and Control (CDC) held jointly with 2009 28th Chinese Control Conference (pp. 3370-3375). IEEE.



Explanation of the files:

"test_subspace.py" - is a driver code for testing the subspace identification method, you should start from here. 

"discretization_test.py" - is used to test the backward Euler method for discretizing a dynamical system

"functionsSID.py" - contains the functions used to estimate the model, the functions from this file are imported in "test_subspace.py" 

All the question and comments should be addressed to me:

Aleksandar Haber 
aleksandar.haber@gmail.com

Also, if you are interested in improving the codes, please let me know. This implementation is not optimized for large-scale problems. 
