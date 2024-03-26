function [PA, PL, X, C, S] = gGKB_gsvd(A, L, k, tol, reorth)
% Compute the first k largest/smallest generalized singular values and right 
% generalized singular vectors of {A,L} by gGKB_gsvd
%
% Reference: [1]. Haibo Li,  Generalizing the Golub-Kahan bidiagonalization for large-scale GSVD computation.
%
% Haibo Li, School of Mathematics and Statistics, The University of Melbourne
% 24, March, 2024.

[m,n] = size(A); 
[p,~] = size(L); 
if (n ~= size(L,2))
  error('The two matrices must have the same column numbers');
end

PA = zeros(m,k);  % save the first k largest left generalized singular vectors
PL = zeros(p,k);  % save the first k largest left generalized singular vectors
X = zeros(n,k);  % save the first k largest right generalized singular vectors
C = zeros(k,1);  % save the first k largest generalized singular values of c_i
S = zeros(k,1);  % save the first k largest generalized singular values of s_i

rng(2024);  % random seed
b1 = randn(m, 1);
b2 = randn(p, 1);

[UA, ZA, ~, BA, ~] = gGKB1_A(A, L, b1, k+1, tol, reorth);
[UL, ~, ~, BL, ~] = gGKB1_L(A, L, b2, k+1, tol, reorth);


BkA = BA(1:k+1,1:k);    
[P1,C1,W1] = svd(BkA, "econ");
C = diag(C1);
S = sqrt(1-C.^2);
PA = UA(:,1:k+1) * P1;   % mx(k+1) x (k+1)xk --> mxk
X  = ZA(:,1:k) * W1;

BkL = BL(1:k+1,1:k);  
[P2,~,~] = svd(BkL, "econ");
% S = diag(S2);
PL = UL(:,1:k+1) * P2;   % px(k+1) x (k+1)xk --> pxk

end
