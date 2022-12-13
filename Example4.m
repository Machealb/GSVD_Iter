%% Accuracy of generalized singular values of JBD-GSVD
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 05, 2022.

clear;
directory = pwd;
path(directory, path)
path([directory, '/Matrices'], path)
path([directory, '/regu'], path)

%% test matrix
N = 800;
c = zeros(N,0);
c(1) = 0.99;   c(2) = 0.98;  c(3) = 0.97;   
c(4:N-3) = linspace(0.96, 0.04, N-6);
c(N-2) = 0.03; c(N-1) = 0.02;  c(N) = 0.01;
s= sqrt(1 - c.*c);
C = diag(c);  S = diag(s);
D = gallery('orthog', N, 2);   
RR = diag(linspace(1,10,N));
invR = diag(1./linspace(1,10,N));
A = C*D'*RR;  L1 = S*D'*RR;
X = invR*D; 

%%
L=1*L1;
[m,n]=size(A);  p=size(L,1);
[Q,R]=qr1([A;L]);  kap = cond(R);
Qa=Q(1:m,:);  Ql=Q(m+1:m+p,:);
k0 = 80;
H1 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at top)
H2 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at bottom)

%--------k0-step JBD algorithm-----------------
tol = 1e-10;
tol2 = 1e-12;
reorth=1;
rng(2022);  % random seed
t = randn(n, 1);
[B, Bh, U, Uh, V_til]=JBD(A,L,t,k0+1,tol,reorth);
[B2, B2h, U2, U2h, V2_til]=JBD(A,L,t,k0+1,tol2,reorth);
P=eye(k0);
for i=1:k0 
    P(i+1,i+1)=-1*P(i,i);
end
V_bar = V_til(:,1:k0+1)*P;
V = Q'*V_til(:,1:k0+1);
Vh = V(:,1:k0+1)*P;
B_bar = Bh(1:k0+1,1:k0+1)*P;
V2_bar = V2_til(:,1:k0+1)*P;
V2 = Q'*V2_til(:,1:k0+1);
V2h = V2(:,1:k0+1)*P;
B2_bar = B2h(1:k0+1,1:k0+1)*P;

%%
er1 = zeros(k0,1);   % error of aproximate c_1 by Bk, tol1
er2 = zeros(k0,1);   % error of aproximate c_1 by Bk, tol2
er3 = zeros(k0,1);   % error of aproximate s_n by Bh_k, tol1
er4 = zeros(k0,1);   % error of aproximate s_n by Bh_k, tol2
bnd = zeros(k0,1);    % upper bound of singular vector error

for k = 1:k0
    fprintf('Running JBD algorithm: the %d-th step ===================\n', k);
    B_k = B(1:k,1:k);    
    [P1,C1,~] = svd(B_k);
    B2_k = B2(1:k,1:k);    
    [P2,C2,~] = svd(B2_k);
    er1(k) = abs(C1(1,1)-c(1));
    er2(k) = abs(C2(1,1)-c(1));
    
    Bh_k = Bh(1:k,1:k);
    [Q1,S1,~] = svd(Bh_k); 
    B2h_k = B2h(1:k,1:k);
    [Q2,S2,~] = svd(B2h_k); 
    er3(k) = abs(S1(k,k)-s(1));
    er4(k) = abs(S2(k,k)-s(1));
end

%------ plot -----------
lw = 2; l = 1:k0;
figure; 
semilogy(l,er1,'->','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.5);
hold on;
semilogy(l,er2,'-s','Color',[0.6350 0.0780 0.1840],'MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',1.5);
hold on;
semilogy(l,kap*tol*ones(k0,1),'-b','LineWidth',2);
hold on;
semilogy(l,kap*tol2*ones(k0,1),'-m','LineWidth',2);
handle = legend('$|c_{1}^{(k)}-c_{1}|, \tau_{1}$','$|c_{1}^{(k)}-c_{1}|, \tau_{2}$',...
    '$\kappa(C)\tau_{1}$','$\kappa(C)\tau_{2}$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
set(gca, 'MinorGridAlpha', 0.02);


figure; 
semilogy(l,er3,'-o','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','[0.3010 0.7450 0.9330]','LineWidth',1.5);
hold on;
semilogy(l,er4,'-^','Color','[0.8500 0.3250 0.0980]','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','[0.8500 0.3250 0.0980]','LineWidth',1.5);
hold on;
semilogy(l,kap*tol*ones(k0,1),'-b','LineWidth',2);
hold on;
semilogy(l,kap*tol2*ones(k0,1),'-m','LineWidth',2);
handle = legend('$|s_{k}^{(k)}-s_{1}|, \tau_{1}$','$|s_{k}^{(k)}-s_{1}|, \tau_{2}$',...
    '$\kappa(C)\tau_{1}$','$\kappa(C)\tau_{2}$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
set(gca, 'MinorGridAlpha', 0.02);
