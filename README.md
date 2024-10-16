# GSVD_Iter

* GSVD_Iter is a MATLAB code library for computing a few Generalized Singular Value Decomposition (GSVD) components of large-scale matrix pairs.

For a matrix pair $\{A,B\}$, if m >= n >= p then the GSVD has the form:

   A = [ U,  0 ] * [ diag(sigma),      0    ] * inv(X)
   
   B = [0, V] *  [ 0,       eye(n-p) \\ \\
                  diag(mu),        0    ] * inv(X)
                    
where

   U  is  m-by-n ,    sigma  is  p-by-1
   
   V  is  p-by-p ,    mu     is  p-by-1
   
   X  is  n-by-n .

Otherwise the GSVD has a more complicated form (see the following reference for more details).
 
Reference: C. F. Van Loan, "Computing the CS and the generalized singular value decomposition", Numerische Mathematik, 46 (1985), pp. 479-491. 

<img src="figs/ritz_A2.png" width="400" />  <img src="figs/ritz_L2.png" width="400" /> 
<img src="figs/conv2_l.png" width="400" /> <img src="figs/conv2_s.png" width="400" />


## JBD_GSVD. 
 
Compute a partial GSVD iteratively using the joint bidigonalization of {A,B}. Here are some research papers related to this method.

1. Zhongxiao Jia, Haibo Li. "[The joint bidiagonalization method for large GSVD computations in finite precision. SIAM Journal on Matrix Analysis and Applications, 44(1), 382–407](https://doi.org/10.1137/22M1483608)."
2. Zhongxiao Jia, Haibo Li. "[The joint bidiagonalization process with partial reorthogonalization. Numerical Algorithm, 88, 965-992](https://doi.org/10.1007/s11075-020-01064-8)."
3. Haibo Li. "[The joint bidiagonalization of a matrix pair with inaccurate inner iterations. SIAM Journal on Matrix Analysis and Applications, 45(1), 232–259](https://doi.org/10.1137/22M1541083)."


## gGKB_GSVD. 
 
Compute a partial GSVD iteratively using the generalized Golub-kahan bidiagonalization. Here are some research papers related to this method.

1. Haibo Li. "[ Characterizing GSVD by singular value expansion of linear operators and its computation](https://arxiv.org/abs/2404.00655)."


## Submit an issue
You are welcome to submit an issue for any questions related to InverProb_IterSolver. 

## License
If you use this code in any future publications, please cite like this:

Zhongxiao Jia, Haibo Li. "[The joint bidiagonalization method for large GSVD computations in finite precision. SIAM Journal on Matrix Analysis and Applications, 44(1), 382–407.](https://doi.org/10.1137/22M1483608)."



