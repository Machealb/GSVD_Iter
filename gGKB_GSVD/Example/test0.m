% verify the basis matrix form recurrences of Q_A and Q_L
clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2023);

%%------------- test matrices---------------
% A = mmread('well1850.mtx');
% L1 = mmread('illc1850.mtx');

% A = mmread('dw2048.mtx');
% L1 = mmread('rdb2048.mtx');

% A = mmread('swang1.mtx');
% [~,n]=size(A);
% L1=get_l(n,1);

% N = 10000;
% c = N:-1:1;  c = c/(2*N);
% s= sqrt(1 - c.*c);
% C = diag(c);  S = diag(s);  
% kapp = 50;
% RR = diag(linspace(1,kapp,N));
% invR = diag(1./linspace(1,kapp,N));
% AA = C*RR;  LL = S*RR;
% A = sparse(AA);  L1 = sparse(LL);
% X = invR; 

N = 500;
c = zeros(N,1);
c(1) = 0.99;   c(2) = 0.98;  c(3) = 0.97;   
c(4:N-3) = linspace(0.96, 0.04, N-6);
c(N-2) = 0.03; c(N-1) = 0.02;  c(N) = 0.01;
% N = nn;
% c = sort(diag(C11(1:nn,1:nn)),'descend');
s= sqrt(1 - c.*c);
C = diag(c);  S = diag(s);
D = gallery('orthog', N, 2);   
RR = diag(linspace(1,10,N));
invR = diag(1./linspace(1,10,N));
A = C*D'*RR;  L1 = S*D'*RR;
X = invR*D; 
m = N; n = N; p = N;

%%---------JBD process------------------------------------------------
L = 1*L1;
[m,n] = size(A);  
p = size(L,1);
M = A'*A + L'*L;

type = 'posi'
k0  = 100;
tol = 0;
reorth = 0;
b = randn(m, 1);
[U, Z, Zb, B, bbeta] = gGKB2_A(A, L, b, k0+1, tol, reorth, type);
% b = randn(p, 1);
% [U, Z, Zb, B, bbeta] = gGKB2_L(A, L, b, k0+1, tol, reorth);

Minv = inv(M);

s1 = zeros(k0, 1);  %
s2 = zeros(k0, 1);  % 
w1 = zeros(k0, 1);  % 
w2 = zeros(k0, 1);  %
bnd1 = zeros(k0,1);
bnd2 = zeros(k0,1);
for k=1:k0
    er1 = A*Z(:,1:k) - U(:,1:k+1)*B(1:k+1,1:k);
    s1(k) = norm(er1);
    Ek = eye(k+1);  ek = Ek(:,k+1);
    er2 = Minv*A'*U(:,1:k+1) - Z(:,1:k)*B(1:k+1,1:k)' - B(k+1,k+1)*Z(:,k+1)*ek';
    s2(k) = norm(er2);

    % er1 = L*Z(:,1:k) - U(:,1:k+1)*B(1:k+1,1:k);
    % s1(k) = norm(er1);
    % Ek = eye(k+1);  ek = Ek(:,k+1);
    % er2 = Minv*L'*U(:,1:k+1) - Z(:,1:k)*B(1:k+1,1:k)' - B(k+1,k+1)*Z(:,k+1)*ek';
    % s2(k) = norm(er2);

    E1 = eye(k+1)-U(:,1:k+1)'*U(:,1:k+1);    w1(k) = norm(E1);
    E2 = eye(k)-Z(:,1:k)'*M*Z(:,1:k);        w2(k) = norm(E2);
end


%%-----------plot--------------------------------------
lw = 1.5; l = 1:k0;

figure; 
semilogy(l,s1,'-b','LineWidth',lw);
hold on;
semilogy(l,s2,'-m','LineWidth',lw);
legend('error1','error2','Location','southeast');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
title('Error of two relations');

figure; 
semilogy(l,w1,'-b','LineWidth',lw);
hold on;
semilogy(l,w2,'-m','LineWidth',lw);
legend('orth1','orth2','Location','southeast');
xlabel('Iteration','Fontsize',15);
ylabel('value','Fontsize',15);
title('Orthogonality');

