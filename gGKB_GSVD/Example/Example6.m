%% Comparison of partial GSVD algorithms

clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2023);

%test matrix
N = 100000;
c = zeros(N,1);
c(1) = 1.0;   c(2) = 0.99;  
c(3:N-2) = linspace(0.98, 0.03, N-4);
c(N-1) = 0.02;  c(N) = 0.01;
s = sqrt(1 - c.^2);
C = spdiags(c(:),0:0,N,N);  
S = spdiags(s(:),0:0,N,N);  
dd = linspace(1,50,N);
RR = spdiags(dd(:),0:0,N,N);
dd1 = 1.0 ./ dd;
invR = spdiags(dd1(:),0:0,N,N);
A1 = C*RR;  L1 = S*RR;
X = invR; 
H1 = speye(N);    H2 = speye(N);    % generalized left singular vectors
m = N; n = N; p = N;


% ----------- basic setting ---------------
A = A1;    L = L1;
k0 = 600;
reorth = 1;  
type = 'posi';
tol = 1e-10;


%%--------- gGKB-GSVD ---------------------
b = randn(m, 1);
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
    x1 = X * [1;zeros(N-1,1)];
    p1 = [1;zeros(N-1,1)];
    ang1_A(k) = angle(p1,p_A1);
    ang1_X(k) = angle(x1,x_k1);

    c_k2 = C1(k,k); 
    w_k2 = W1(:,k);   
    x_k2 = Z(:,1:k) * w_k2;
    p_A2 = U(:,1:k+1) * P1(:,k);
    er2(k)  = abs(c(n)-c_k2);
    xn = X * [zeros(N-1,1);1];
    pn = [zeros(N-1,1);1];
    ang2_A(k) = angle(pn,p_A2);
    ang2_X(k) = angle(xn,x_k2);
end


%%--------- JBD-GSVD ---------------------
t = ones(n, 1);
[B1, Bh1, U1, Uh1, V_til1]=JBD(A,L,t,k0+1,tol,reorth);
P=eye(k0);
for i=1:k0 
    P(i+1,i+1)=-1*P(i,i);
end
V_bar1 = V_til1(:,1:k0+1)*P;
B_bar1 = Bh1(1:k0+1,1:k0+1)*P;

er3 = zeros(k0,1);        % error of c_1
ang3_A = zeros(k0,1);     % angle error of P_A1
ang3_X = zeros(k0,1);     % angle error of X_1
er4 = zeros(k0,1);        % error of c_n
ang4_A = zeros(k0,1);     % angle error of P_An
ang4_X = zeros(k0,1);     % angle error of X_n


for k = 1:k0
    fprintf('Running JBD algorithm: the %d-th step ===================\n', k);
    B_k1 = B1(1:k,1:k);    
    [P1,C1,W1] = svd(B_k1);

    c_k1 = C1(1,1); 
    w_k1 = W1(:,1);
    p_A1 = U1(:,1:k) * P1(:,1);  
    x_k1 = lsqr(@(x,tflag)afun(x,A,L,tflag),V_til1(:,1:k)*w_k1,1e-8,n);
    er3(k)  = abs(c(1)-c_k1);
    % fprintf('Norm of p_A1: %f\n', norm(p_A1));
    x1 = X * [1;zeros(N-1,1)];
    p1 = [1;zeros(N-1,1)]; 
    ang3_X(k) = angle(x1,x_k1);
    ang3_A(k) = angle(p1,p_A1);
    

    c_k2 = C1(k,k); 
    w_k2 = W1(:,k);  
    p_A2 = U1(:,1:k) * P1(:,k);
    x_k2 = lsqr(@(x,tflag)afun(x,A,L,tflag),V_til1(:,1:k)*w_k2,1e-8,n);
    er4(k)  = abs(c(n)-c_k2);
    xn = X * [zeros(N-1,1);1];
    pn = [zeros(N-1,1);1];
    ang4_X(k) = angle(xn,x_k2);
    ang4_A(k) = angle(pn,p_A2);

    
end


%------ plot ---------------------------
l = 1:k0;

% figure('units','normalized','position',[0.1,0.1,1.0,0.8]);

figure;
semilogy(l(1:100),er1(1:100),'-v','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l(1:100),er3(1:100),'-o','Color','[0.8500 0.3250 0.0980]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.8500 0.3250 0.0980]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $c_{1}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l(1:100),ang1_A(1:100),'-v','Color','[0.4660 0.6740 0.1880]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.4660 0.6740 0.1880]','LineWidth',1.5);
hold on;
semilogy(l(1:100),ang3_A(1:100),'-o','Color','[0.4940 0.1840 0.5560]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.4940 0.1840 0.5560]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $p_{A,1}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l(1:100),ang1_X(1:100),'-v','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l(1:100),ang3_X(1:100),'-o','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:9:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.3010 0.7450 0.9330]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $x_{1}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);


figure;
semilogy(l(1:600),er2(1:600),'-v','Color','b','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l(1:600),er4(1:600),'-o','Color','[0.8500 0.3250 0.0980]','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.8500 0.3250 0.0980]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $c_{n}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l(1:600),ang2_A(1:600),'-v','Color','[0.4660 0.6740 0.1880]','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.4660 0.6740 0.1880]','LineWidth',1.5);
hold on;
semilogy(l(1:600),ang4_A(1:600),'-o','Color','[0.4940 0.1840 0.5560]','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.4940 0.1840 0.5560]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $p_{A,n}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure;
semilogy(l(1:600),ang2_X(1:600),'-v','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l(1:600),ang4_X(1:600),'-o','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:19:k0,...
'MarkerSize',5,'MarkerFaceColor','[0.3010 0.7450 0.9330]','LineWidth',1.5);
legend('gGKB\_GSVD', 'JBD\_GSVD','Fontsize',18,'interpreter','latex');
xlabel('Iteration','Fontsize',18);
ylabel('Error','Fontsize',18);
title('Convergence for $x_{n}$','Fontsize',19,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);
