function [B, Bh, U, Uh, V_til] = JBD(A,L,t,k,tol,reorth)
%  Joint bidiagonalization reduction for matrix pair {A,B}, Zha's version
%  Starting vector is t 
%  Reduce A to upper bidiagonal matrix B of kxk
%  Reduce L to upper bidiagonal matrix Bh of kxk
%  where U, Uh and V consist of the left and right Lanczos vectors.
%
%  QR factorization of [A;L] is avoided by iteratively solving a least squares 
%  problems with matrix [A;L], using LSQR solver with tolerance tol
%
%  Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization,
%    reorth = 1 : reorthogonalization of v_til by means of MGS,
%    reorth = 2 : double reorthogonalization of v_til
%
% Reference: H. Zha, "Computing the generalized singular 
% values/vectors of large sparse or structured matrix pairs",
% Numerische Mathematik, 72 (1996), pp. 391–417.
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 03, 2022.

% Initialization.
[m,n] = size(A); p=size(L,1);
if (n ~= size(L,2))
  error('The two matrices must have the same column numbers')
end
C = [A;L];  r = C*t;
v_til = r/norm(r);  V_til(:,1) = v_til;
alpha = norm(v_til(1:m));  B(1,1) = alpha;
u = v_til(1:m)/alpha;  U(:,1) = u;
alpha_h = norm(v_til(m+1:m+p)); Bh(1,1) = alpha_h;
u_h = v_til(m+1:m+p)/alpha_h;  Uh(:,1) = u_h;

% The k-step iteration
for i=1:k
    u_til=[u; zeros(p,1)];
    z = lsqr(@(x,tflag)afun(x,A,L,tflag),u_til,tol,n);  
    r = C*z - alpha*v_til;

    % Reorthogonalization of v_til
    if reorth == 0
        r=r;
    elseif reorth == 1
        for j=1:i
            r = r - V_til(:,j)*V_til(:,j)'*r;
        end
    elseif reorth == 2
        for j=1:i  r = r - (V_til(:,j)'*r)*V_til(:,j);  end
        for j=1:i  r = r - (V_til(:,j)'*r)*V_til(:,j);  end
    end
    beta = norm(r);  B(i,i+1) = beta;
    v_til = r/beta;  V_til(:,i+1) = v_til;
    beta_h = alpha*beta/alpha_h;  Bh(i,i+1) = beta_h;

    % Compute u and u_h
    s = v_til(1:m) - beta*u;
    alpha = norm(s); B(i+1,i+1) = alpha;
    u = s/alpha;  U(:,i+1) = u;  
    if mod(i,2)==0
        vv = v_til(m+1:m+p);
    else
        vv = -v_til(m+1:m+p);
    end
    q = vv - beta_h*u_h;
    alpha_h = norm(q); Bh(i+1,i+1) = alpha_h;
    u_h = q/alpha_h;  Uh(:,i+1) = u_h;
end

end



