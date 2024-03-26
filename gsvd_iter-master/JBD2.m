function [B, Bbar, U, U_hat, V_til, bbeta]=JBD2(A,L,b,k,tol,reorth)
%  Joint bidiagonalization reduction for matrix pair {A,B}， Kilmer's version
%  with starting vector b, and bbeta=||b||
%  Reduce A to lower bidiagonal matrix B of (k+1)xk
%  Reduce L to upper bidiagonal matrix Bh of kxk
%  where U, U_hat and V_til consist of the left and right Lanczos vectors.
%
%  QR factorization of [A;L] is avoided by iteratively solving a least squares 
%  problems with matrix [A;L], using LSQR solver with tolerance tol
%
%  Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization,
%    reorth = 1 : reorthogonalization by means of MGS,
%    reorth = 2 : double reorthogonalization
% 
% Reference: M. E. Kilmer, P. C. Hansen, and M. I. Espanol, "A projection-based 
% approach to general-form Tikhonov regularization",
% SIAM J. Sci. Comput., 29 (2007), pp. 315–330.
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 03, 2022.

% Initialization.
[m,n] = size(A); p=size(L,1);
if (n ~= size(L,2))
  error('The two matrices must have the same column numbers')
end
beta=norm(b);  bbeta=beta;
u=b/beta;  U(:,1)=u;

utilde=[u; zeros(p,1)];
x = lsqr(@(z,tflag)afun(z,A,L,tflag),utilde,tol,n);
ss = A*x; tt = L*x;
v = [ss;tt];  alpha=norm(v);
v = v/alpha;  B(1,1) = alpha;  V_til(:,1) = v;

uhat = v(m+1:m+p);  alphahat = norm(uhat);
uhat = uhat/alphahat;  Bbar(1,1) = alphahat;  U_hat(:,1) = uhat;

if (reorth == 0)
    u = v(1:m) - alpha * u;
elseif (reorth == 1)
    u = v(1:m) - alpha * u;
    u = u - U * (U' * u);
elseif (reorth == 2)
    u = v(1:m) - alpha * u; 
    u = u - U * (U' * u);
    u = u - U * (U' * u);
end
beta = norm(u);  u = u/beta;  B(2,1) = beta;  U(:,2) = u;

% the k-step iteration
for i = 2:k
    utilde = [U(:,i); zeros(p,1)];
    x = lsqr(@(z,tflag)afun(z,A,L,tflag), utilde,tol,n);
    ss = A*x; tt = L*x;
    Qu = [ss;tt];
    if (reorth == 0)
        v = Qu - B(i, i-1)*V_til(:,i-1);
    elseif(reorth == 1)
        v = Qu - B(i, i-1)*V_til(:,i-1);
        for j=1:i-1, v = v - (V_til(:,j)'*v)*V_til(:,j); end
    elseif (reorth == 2)
        v = Qu - B(i, i-1)*V_til(:,i-1);
        for j=1:i-1, v = v - (V_til(:,j)'*v)*V_til(:,j); end
        for j=1:i-1, v = v - (V_til(:,j)'*v)*V_til(:,j); end
    end
    alpha = norm(v);
    v = v/alpha;
    B(i,i) = alpha;
    V_til(:,i) = v;
    
    betahat=(alpha*B(i,i-1))/alphahat;
    if(mod(i,2)==0)
        Bbar(i-1,i) = -betahat;
    else
        Bbar(i-1,i) = betahat;
    end
    
    if(mod(i,2)==0)
        vv = -v(m+1:m+p);
    else
        vv = v(m+1:m+p);
    end
    if (reorth == 0)
        uhat = vv - betahat * U_hat(:,i-1);
    elseif (reorth == 1)
        uhat = vv - betahat * U_hat(:,i-1);
        for j=1:i-1, uhat = uhat - (U_hat(:,j)'*uhat)*U_hat(:,j); end
    elseif (reorth == 2)
        uhat = vv - betahat * U_hat(:,i-1);
        for j=1:i-1, uhat = uhat - (U_hat(:,j)'*uhat)*U_hat(:,j); end
        for j=1:i-1, uhat = uhat - (U_hat(:,j)'*uhat)*U_hat(:,j); end
    end
    alphahat = norm(uhat);
    if(mod(i,2)==0)
        Bbar(i,i) = -alphahat;
    else
        Bbar(i,i) = alphahat;
    end
    uhat = uhat/alphahat;
    U_hat(:,i) = uhat;
    
    if (reorth == 0)
        u = v(1:m) - alpha * u;
    elseif (reorth == 1)
        u = v(1:m) - alpha * u;
        for j=1:i, u = u - (U(:,j)'*u)*U(:,j); end
    elseif (reorth == 2)
        u = v(1:m) - alpha * u;
        for j=1:i, u = u - (U(:,j)'*u)*U(:,j); end
        for j=1:i, u = u - (U(:,j)'*u)*U(:,j); end
    end
    beta = norm(u);
    u = u/beta;
    B(i+1,i) = beta;
    U(:,i+1) = u;
end

end


