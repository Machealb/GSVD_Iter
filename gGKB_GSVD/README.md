# gGKB_GSVD
Using generalizing Golub-Kahan bidiagonalization to compute a part of nontrivial GSVD components.

Plase read the following guideline to run the codes and test the algorithms.

## test recursive relations and orthogonality
test0.m

## 1. Convergence of c_i and s_i, tol = 0, reorth = 1 or 0

gGKB of A and gGKB of L, respectively; {A, L} is self-constructed

Example1.m;  Example2.m.

## 2. Convergence of GSVD components, residual norm and upper bound, reorth = 1, tol=0

Show the gGKB of A，for largest and smallest GSVD components.
{A, L} are rectangles, compute its GSVD by matlab's gsvd function as baseline.
{well1850, L1}.

Example3.m.


## 3. M is positive semidefinite, convergence and accuracy of GSVD components

Show the gGKB of A，for largest and smallest GSVD components;
{A, L} is self-constructed.

Example4.m.

## 4. Final accuracy influenced by the value of tol

Set a relative larger-scale {A, L}; use tol=1e-10 and 1e-8.

Example5.m

