# Tests for gGKB_GSVD

## test recursive relations and orthogonality
test0.m


## Show1. ghost convergence of c_i and s_i, smallest and largest, tol = 0, reorth = 1 or 0

gGKB of A and gGKB of L, respectively:
{A, L} is self-constructed

Example1.m:
ritz_A1   Ritz value without reorthgonalization
ritz_A2   Ritz value with reorthgonalization
erA_c1   error for first 3 c_i
erA_c2    error for last  3 c_i

Example2.m:
ritz_L1   Ritz value without reorthgonalization
ritz_L2   Ritz value with reorthgonalization
erL_s1    error for first 3 s_i
erL_s2    error for last  3 s_i


## Show2. Convergence of GSVD components, residual and upper bound, full reorthogonaliztion, tol=0

Only for gGKB of A， largest and smallest GSVD components.
{A, L} are rectangles, compute its GSVD by matlab's gsvd function.
{well1850, L1}

Example3.m:
conv1_l   convergence of largest  GSVD components
conv1_s   convergence of smallest GSVD components


## Show3. M is positive semidefinite, convergence and accuracy of GSVD components

Only for gGKB of A, largest and smallest GSVD components.
{A, L} is self-constructed

Example4.m:
conv2_l   convergence of largest  GSVD components
conv2_s   convergence of smallest GSVD components

## Show4. Final accuracy influenced by the value of tol

Choose a relative larger-scale {A, L}, compute its GSVD by matlab's gsvd function.
{dw2048, rdb2048} 

set reorth = 1

Example5.m: 
conv3_l   convergence of largest  GSVD components, tol=1e-10
conv3_s   convergence of smallest GSVD components, tol=1e-10
conv4_l   convergence of largest  GSVD components, tol=1e-8
conv4_s   convergence of smallest GSVD components, tol=1e-8




