function [X,C,S] = jbdgsvd(A,B,k,tol,K)
% Compute the first k largest generalized singular values and right 
% generalized singular vectors of {A,B} by rJBD
%
% The QR factorization of [A;B] is avoided by iteratively solving a least squares 
% problems with matrix [A;B], using LSQR solver with tolerance tol
% 
% The maximun iteration is K.


[~,n] = size(A); 
if (n ~= size(B,2))
  error('The two matrices must have the same column numbers');
end

X = zeros(n,1);  % save the first k largest right generalized singular vectors
C = zeros(k,1);  % save the first k largest generalized singular values of c_i
S = zeros(k,1);  % save the first k largest generalized singular values of s_i

rng(2022);  % random seed
t = randn(n, 1);
tol2 = 10 * tol; % stopping tolerance for get x_i
[B, Bh, ~, ~, V_til] = JBD(A,B,t,K,tol,1);
for l = k:K
    B_l = B(1:l,1:l);    
    [~,C1,W1] = svd(B_l);
    C = C1(1:k,1:k);
    S = sqrt(1-C.^2);
    for i = 1:k
        X(:,i) = lsqr(@(x,tflag)afun(x,A,L,tflag),V_til(:,1:k)*W1(:,i),tol2,n);
    end
end

end