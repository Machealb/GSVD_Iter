%% Accuracy of GSVD components by gGKB_GSVD,
% Ms = v is solved inaccurately.

clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2023);

% ----- test matrix -----
N = 5000;
c = zeros(N,1);
c(1) = 0.99;   c(2) = 0.97;  
c(3:N-2) = linspace(0.95, 0.15, N-4);
c(N-1) = 0.1;  c(N) = 0.05;
s = sqrt(1 - c.^2);
C = diag(c);  S = diag(s);
D = gallery('orthog', N, 2);   
RR = diag(linspace(1,10,N));
invR = diag(1./linspace(1,10,N));
A1 = C*D'*RR;  L1 = S*D'*RR;
X = invR*D; 
H1 = eye(N);    H2 = eye(N);    % generalized left singular vectors
m = N; n = N; p = N;


%%--------- gGKB set up ---------------------
A = A1;     L = L1;
AA = A'*A;  LL = L'*L; 
kap = sqrt(norm(AA+LL));

% k0-step gGKB algorithm
type = 'posi';
tol = 1e-10;
b = randn(m, 1);
k0 = 150;
reorth = 1; 
[U, Z, Zb, B, bbeta] = gGKB2_A(A, L, b, k0+1, tol, reorth, type);
 

er1 = zeros(k0,1);        % error of c_1
ang1_A = zeros(k0,1);     % angle error of P_A1
ang1_X = zeros(k0,1);     % angle error of X_1
er2 = zeros(k0,1);        % error of c_n
ang2_A = zeros(k0,1);     % angle error of P_An
ang2_X = zeros(k0,1);     % angle error of X_n

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
    ang1_X(k) = angle(X(:,1),x_k1);

    c_k2 = C1(k,k); 
    w_k2 = W1(:,k);   
    x_k2 = Z(:,1:k) * w_k2;
    p_A2 = U(:,1:k+1) * P1(:,k);
    er2(k)  = abs(c(n)-c_k2);
    ang2_A(k) = angle(H1(:,n),p_A2);
    ang2_X(k) = angle(X(:,n),x_k2);
end


%------ plot ---------------------------
l = 1:k0;

figure;
semilogy(l(1:100),er1(1:100),'-*','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l(1:100),ang1_A(1:100),'-o','Color','m','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','m','LineWidth',1.5);
hold on;
semilogy(l(1:100),ang1_X(1:100),'->','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','g','LineWidth',1.5);
handle = legend('$|c_{1}-\theta_{1}^{(k)}|$','$\sin\angle(p_{A,1},\bar{p}_{A,1}^{(k)})$',...
    '$\sin\angle(x_{1},\bar{x}_{1}^{(k)})$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The 1-st GSVD components','Fontsize',16,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l,er2,'-*','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l,ang2_A,'-o','Color','m','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','m','LineWidth',1.5);
hold on;
semilogy(l,ang2_X,'->','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','g','LineWidth',1.5);
handle = legend('$|c_{n}-\theta_{k}^{(k)}|$','$\sin\angle(p_{A,n},\bar{p}_{A,k}^{(k)})$',...
    '$\sin\angle(x_{n},\bar{x}_{k}^{(k)})$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The n-th GSVD components','Fontsize',16,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);



