% show the error f_{l+1} of recrrence of [\beta_l, \alpha_{l+1}]
% show the relation of B and Bh by the computation error
% show the error of reduction of Ql
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 04, 2022.

clear;
directory = pwd;
path(directory, path)
path([directory, '/Matrices'], path)

%%------------- test matrices---------------
% A = mmread('illc1850.mtx');
% L1 = mmread('well1850.mtx');

A = mmread('swang1.mtx');
[~,n]=size(A);
L1=get_l(n,1);

% 
% N = 600;
% c = zeros(N,0);
% c(1) = 0.99;   c(2) = 0.98;  c(3) = 0.97;   
% c(4:N-3) = linspace(0.96, 0.04, N-6);
% c(N-2) = 0.03; c(N-1) = 0.02;  c(N) = 0.01;
% % N = nn;
% % c = sort(diag(C11(1:nn,1:nn)),'descend');
% s= sqrt(1 - c.*c);
% C = diag(c);  S = diag(s);
% D = gallery('orthog', N, 2);   
% RR = diag(linspace(1,10,N));
% invR = diag(1./linspace(1,100,N));
% A = C*D'*RR;  L1 = S*D'*RR;
% %X = invR*D; 
% m = N; n = N; p = N;

%%---------------------------------------------------------
L = 0.1*L1;
[m,n]=size(A);  p=size(L,1);
C = [A;L];  
[Q,R]=qr1(C);  kap = cond(R);
Qa=Q(1:m,:);  Ql=Q(m+1:m+p,:);
% QQ = [C; S]; Q = sparse(QQ);  kap = kapp;
% Qa = sparse(C);  Ql = sparse(S);

k0=60;
tol = 1e-8;
reorth=1;
rng(2022);  % random seed
%t = randi(10,n,1);
t = randn(n, 1);

[B, Bh, U, Uh, V_til]=JBD(A,L,t,k0+1,tol,reorth);
P=eye(k0+1);
for i=1:k0 
    P(i+1,i+1)=-1*P(i,i);
end
V_bar = V_til(:,1:k0+1)*P;
V = Q'*V_til(:,1:k0+1);
Vh = V(:,1:k0+1)*P;
B_bar = Bh(1:k0+1,1:k0+1)*P;
s1 = zeros(k0, 1);  % error of recurrence of [\beta_l,\alpha_{l+1}]
bnd1 = zeros(k0,1);  % upper bound of s1
s2 = zeros(k0, 1);  % error of relation of B and Bh
bnd2 = zeros(k0,1);  % upper bound of s2
s3 = zeros(k0, 1);  % error of reduction of Ql
bnd3 = zeros(k0,1);  % upper bound of s3

ey1 = eye(n);
p_1 = [-ey1(:,1); U(:,1)];
P_l = eye(m+n) - p_1*p_1';
for l = 1:k0-1
    %fprintf('Iteration %d\n', l);
    p_l = [-ey1(:,l+1); U(:,l+1)];
    P_l = P_l * (eye(m+n) - p_l*p_l');
    El = eye(l);  el = B(l,l+1)*El(:,l);
    Es = eye(m+n-l);  es = B(l+1,l+1)*Es(:,1);
    B_sl = [el; es];  % [\beta_l, \alpha_{l+1}]
    f_l1 = P_l*B_sl - [zeros(n,1); Qa*V(:,l+1)];  % f_{l+1}
    s1(l+1) = norm(f_l1);
end
bnd1= 5*kap*tol*ones(k0,1);

for k = 1:k0
    er1 = eye(k) - B(1:k,1:k)'*B(1:k,1:k) - B_bar(1:k,1:k)'*B_bar(1:k,1:k);
    s2(k) = norm(er1);
    bnd2(k) = 5*kap*tol;
    Ek = eye(k);  ek = Ek(:,k);
    er2 = Ql'*Uh(:,1:k) - Bh(k,k+1)*Vh(:,k+1)*ek';
    s3(k) = norm(Vh(:,k+1)'*er2);  % ||\hat{v}_{k+1}^{T}\hat{G}_k||
    sig = svd(Bh(1:k,1:k));
    bnd3(k) = kap*tol/sig(k);
end

%----------- plot ------------------
lw = 2; l1 = 1:k0;
figure; 
semilogy(l1,s1,'-.r','LineWidth',2);
hold on;
semilogy(l1,s2,'-.g','LineWidth',2);
hold on;
semilogy(l1,bnd1,'-b','LineWidth',2)
handle=legend('$\|f_{l+1}\|$','$\|H_{k}\|$','$5\kappa(C)\tau$','Location','southeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
%set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.02);


% l2 = 1:k0;
% figure; 
% semilogy(l2,s2,'-.b','LineWidth',lw);
% hold on;
% semilogy(l2,bnd2,'-r','LineWidth',lw);
% handle=legend('$\|H_{k+1}\|$','$l\kappa(C)\tau$','Location','southeast');
% set(handle,'interpreter','latex');
% xlabel('Iteration','Fontsize',12);
% ylabel('Error','Fontsize',12);
% grid on;

figure; 
semilogy(l1,s3,'-b','LineWidth',2);
hold on;
semilogy(l1,bnd3,'-r','LineWidth',2)
handle=legend('$\|\hat{v}_{k+1}^{T}G_{k}\|$','$\|\widehat{B}_{k}^{-1}\|\kappa(C)\tau$','Location','southeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
%set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.02);
    