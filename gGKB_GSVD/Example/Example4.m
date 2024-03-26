%% Accuracy of GSVD components by gGKB_GSVD
% M = A'*A + L'*L is  positive semidefinite

clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2023);

%% test matrix
N = 1000;
r0 = 100;
r = N - r0;
c = zeros(N,1);
c(1) = 0.99;   c(2) = 0.98;  
c(3:r-2) = linspace(0.96, 0.06, r-4);
c(r-1) = 0.04;  c(r) = 0.02;
s = zeros(N,1);
cr = c(1:r);
s(1:r) = sqrt(1 - cr.^2);
C = diag(c);  S = diag(s);
D = gallery('orthog', N, 2);   
RR = diag(linspace(1,10,N));
invR = diag(1./linspace(1,10,N));
A1 = C*D'*RR;  L1 = S*D'*RR;
X = invR*D; 
H1 = eye(N);    H2 = eye(N);    % generalized left singular vectors
m = N; n = N; p = N;

A = A1;     L = L1;
AA = A'*A;  LL = L'*L;
M = AA + LL;
X1 = M * pinv(M) * X(:,1:r);


%% gGKB set up
kap = sqrt(norm(M));
type = 'semi';
tol = 0;
b = randn(m, 1);
k0 = 150;
reorth = 1;  
[U, Z, Zb, B, bbeta] = gGKB2_A(A, L, b, k0+1, tol, reorth, type);

er1 = zeros(k0,1);        % error of c_1
er2 = zeros(k0,1);        % error of c_n
ang1_A = zeros(k0,1);     % angle error of P_A1
ang1_X = zeros(k0,1);     % angle error of X_1
ang2_A = zeros(k0,1);     % angle error of P_An
ang2_X = zeros(k0,1);     % angle error of X_n
res1_nm = zeros(k0,1);    % residual norm, for c_1
res1_bn = zeros(k0,1);    % residual norm bound, for c_1
res2_nm = zeros(k0,1);    % residual norm, for c_n
res2_bn = zeros(k0,1);    % residual norm bound, for c_n

for k = 1:k0
    fprintf('Running gGKB_GSVD algorithm: the %d-th step ===================\n', k);
    B_k = B(1:k+1,1:k);    
    [P1,C1,W1] = svd(B_k, "econ");
    c_k1 = C1(1,1);      
    w_k1 = W1(:,1);
    x_k1 = Z(:,1:k) * w_k1;
    p_A1 = U(:,1:k+1) * P1(:,1);
    er1(k)  = abs(c(1)-c_k1);
    ang1_A(k) = angle(H1(:,1),p_A1);
    ang1_X(k) = angle(X1(:,1),x_k1);
    res1= ((1-c_k1^2)*AA-c_k1^2*LL) * x_k1;
    res1_nm(k) = norm(res1) / kap;
    res1_bn(k) = abs(B(k+1,k+1)*B(k+1,k)*w_k1(k))+50*eps;

    c_k2 = C1(k,k); 
    w_k2 = W1(:,k);   
    x_k2 = Z(:,1:k) * w_k2;
    p_A2 = U(:,1:k+1) * P1(:,k);
    er2(k)  = abs(c(r)-c_k2);
    ang2_A(k) = angle(H1(:,r),p_A2);
    ang2_X(k) = angle(X1(:,r),x_k2);
    res2= ((1-c_k2^2)*AA-c_k2^2*LL) * x_k2;
    res2_nm(k) = norm(res2) / kap;
    res2_bn(k) = abs(B(k+1,k+1)*B(k+1,k)*w_k2(k))+10*eps;

end


%------ plot ---------------------------
l = 1:k0;

figure;
semilogy(l,er1,'-v','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l,ang1_A,'-s','Color','r','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','r','LineWidth',1.5);
hold on;
semilogy(l,ang1_X,'->','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.3010 0.7450 0.9330]','LineWidth',1.5);
hold on;
semilogy(l,res1_nm,'-*','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l,res1_bn,'-o','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',5,'LineWidth',1.5);
handle = legend('$|c_{1}-\theta_{1}^{(k)}|$','$\sin\angle(p_{A,1},\bar{p}_{A,1}^{(k)})$',...
    '$\sin\angle(x_{1},\bar{x}_{1}^{(k)})$','$\|r_{1}^{(k)}\|_2/\|M\|_{2}^{1/2}$',...
    'upper bound','Location','northeast');
set(handle,'Fontsize',15,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The 1-st GSVD components','Fontsize',16,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l,er2,'-v','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l,ang2_A,'-D','Color','r','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','r','LineWidth',1.5);
hold on;
semilogy(l,ang2_X,'->','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.3010 0.7450 0.9330]','LineWidth',1.5);
hold on;
semilogy(l,res2_nm,'-*','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l,res2_bn,'-o','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',5,'LineWidth',1.5);
handle = legend('$|c_{r}-\theta_{k}^{(k)}|$','$\sin\angle(p_{A,r},\bar{p}_{A,k}^{(k)})$',...
    '$\sin\angle(x_{r},\bar{x}_{k}^{(k)})$','$\|r_{k}^{(k)}\|_2/\|M\|_{2}^{1/2}$',...
    'upper bound','Location','northeast');
set(handle,'Fontsize',15,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The r-th GSVD components','Fontsize',16,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);
