%% Show the ghost generalized singular values phenomenon by gGKB-GSVD_L algorithm

clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2023);

% ----- test matrix -----
N = 1000;
c = zeros(N,1);
c(1) = 1;   c(2) = 0.95;  c(3) = 0.90;   
c(4:N-3) = linspace(0.88, 0.12, N-6);
c(N-2) = 0.10; c(N-1) = 0.05;  c(N) = 0.01;
s= sqrt(1 - c.*c);
C = diag(c);  S = diag(s);  
kapp = 100;
RR = diag(linspace(1,kapp,N));
invR = diag(1./linspace(1,kapp,N));
AA = C*RR;  LL = S*RR;
A1 = sparse(AA);  L1 = sparse(LL);
X = invR;
m = N; n = N; p = N;

%%
A = A1;     L = L1;
k0 = 150;
H1 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at top)
H2 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at bottom)

%%  k0-step JBD algorithm
type = 'posi';
tol = 0;
reorth=0;
b = randn(m, 1);
[U, Z, Zb, B, bbeta] = gGKB2_L(A, L, b, k0+1, tol, reorth, type);

for k = 1:k0
    fprintf('Running sub_SVD: the %d-th step -------\n', k);
    B_k = B(1:k+1,1:k);
    s_k = svd(B_k);  
    H1(1:k, k) = s_k;  % approximate  gen singular values S (largest)
    H2(k0-k+1:k0, k) = s_k;  % approximate gen singular values S (smallest)

end

%% plotting
figure;
for i = 1: 3
    x_i = i : 1 : k0;
    plot(x_i, H1(i, i:k0), 'x-');
    hold on;
end
scatter(k0*ones(N, 1), s, 'ko');
xlabel('Iteration','Fontsize',15);
ylabel('Ritz value','Fontsize',15);

for i = 1: 3
    x_i = i : 1 : k0;
    plot(x_i, H2(k0-i+1, i:k0), 'x-');
    hold on;
end
scatter(k0*ones(N, 1), s, 'ko');
xlabel('Iteration','Fontsize',15);
ylabel('Ritz  value','Fontsize',15);
title('Approximations to $s_{i}$', 'Fontsize',17,'interpreter','latex');


%% plot relaltive error w.r.t. k
er1 = zeros(k0,1);
er2 = zeros(k0,1);
er3 = zeros(k0,1);
er4 = zeros(k0,1);
er5 = zeros(k0,1);
er6 = zeros(k0,1);

for i = 1:k0 
    er1(i) = abs(s(N)-H1(1,i));   % error of s_N
    er2(i) = abs(s(N-1)-H1(2,i));   % error of s_{N-1}
    er3(i) = abs(s(N-2)-H1(3,i));   % error of s_{N-2}
    
    er4(i) = abs(s(2)-H2(k0,i));  % error of s_1
    er5(i) = abs(s(3)-H2(k0-1,i));  % error of s_2
    er6(i) = abs(s(4)-H2(k0-2,i));  % error of s_3
end

er2(1) = 0;  
er3(1) = 0;  er3(2) = 0;
er5(1) = 0;
er6(1) = 0;  er6(2) = 0;

lw = 2.0; l = 1:k0;

figure; 
semilogy(l,er1,'-d','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l,er2,'-o','Color','b','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l,er3,'-v','Color','m','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','m','LineWidth',1.5);
handle=legend('$|s_{n}-\theta_{k}^{(k)}|$','$|s_{n-1}-\theta_{k-1}^{(k)}|$', '$|s_{n-2}-\theta_{k-2}^{(k)}|$','Location','northeast');
set(handle,'Fontsize',16,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The last three $s_{i}$', 'Fontsize',17,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure; 
semilogy(l,er4,'-d','Color','[0.6350 0.0780 0.1840]','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','[0.6350 0.0780 0.1840]','LineWidth',1.5);
hold on;
semilogy(l,er5,'-o','Color','b','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','b','LineWidth',1.5);
hold on;
semilogy(l,er6,'-v','Color','m','MarkerIndices',1:9:k0,...
    'MarkerSize',5,'MarkerFaceColor','m','LineWidth',1.5);
handle=legend('$|s_{2}-\theta_{1}^{(k)}|$','$|s_{3}-\theta_{2}^{(k)}|$', '$|s_{4}-\theta_{3}^{(k)}|$','Location','northeast');
set(handle,'Fontsize',16,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('The first three $s_{i}$', 'Fontsize',17,'interpreter','latex');
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);



