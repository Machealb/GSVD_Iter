# gsvd_iter
Generalized singular value decomposion (GSVD) computation by iterative methods.

For a matrix pair $\{A,B\}$, if m >= n >= p then the GSVD has the form:
   [ A ] = [ U  0 ]*[ diag(sigma)      0    ]*inv(X)
   [ B ]   [ 0  V ] [      0       eye(n-p) ]
                    [  diag(mu)        0    ]
where
   U  is  m-by-n ,    sigma  is  p-by-1
   V  is  p-by-p ,    mu     is  p-by-1
   X  is  n-by-n .

Otherwise the GSVD has a more complicated form (see manual for details).
 
Reference: C. F. Van Loan, "Computing the CS and the generalized 
singular value decomposition", Numer. Math. 46 (1985), 479-491. 
