%% Accuracy of generalized right vectors of JBD-GSVD, by error norm
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 05, 2022.

clear;
directory = pwd;
path(directory, path)
path([directory, '/Matrices'], path)

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
tol = 1e-12;
reorth=1;
rng(2022);  % random seed
t = randn(n, 1);
[B, Bh, U, Uh, V_til]=JBD(A,L,t,k0+1,tol,reorth);
P=eye(k0);
for i=1:k0 
    P(i+1,i+1)=-1*P(i,i);
end
V_bar = V_til(:,1:k0+1)*P;
V = Q'*V_til(:,1:k0+1);
Vh = V(:,1:k0+1)*P;
B_bar = Bh(1:k0+1,1:k0+1)*P;

%%
tol11 = 1*tol;  % tolerance1 of solving Cx=b
tol12 = 100*tol;  % tolerance2 of solving Cx=b
ang0 = zeros(k0,1);  % norm error of x1 using QR
er0 = zeros(k0,1);   % norm error of x1 using QR
er1 = zeros(k0,1);   % norm error of x1 using tol1
er2 = zeros(k0,1);   % norm error of x1 using tol2
bnd = zeros(k0,1);   % upper bound of singular vector error
for k = 1:k0
    fprintf('Running JBD algorithm: the %d-th step ===================\n', k);
    B_k = B(1:k,1:k);    
    [P1,C1,W1] = svd(B_k);
    % B_hk = B_hat(1:k,1:k);
    % [Q1,S1,W2] = svd(B_hk); 
    % 
    w_k1 = W1(:,1);  
    x_k0 = R \ (Q'*V_til(:,1:k)*w_k1);  % first approximate generalized singular vector
    x_k1 = lsqr(@(x,tflag)afun(x,A,L,tflag),V_til(:,1:k)*w_k1,tol11,n);
    x_k2 = lsqr(@(x,tflag)afun(x,A,L,tflag),V_til(:,1:k)*w_k1,tol12,n);
    ang0(k) = angle(X(:,1),x_k0);   % sin(angle) of right gen-singular vector x_k0
    er0(k) = norm(X(:,1)-x_k0)/1.0; % ||x_1-x_k1||/norm(inv(R)), 
    er1(k) = norm(X(:,1)-x_k1)/1.0;   % right gen-singular vector x_k2
    er2(k) = norm(X(:,1)-x_k2)/1.0;   %  right gen-singular vector x_k2
    bnd(k) = kap*tol/(c(1)-c(2));
end

%------ plot -----------
l = 1:k0;
figure;
semilogy(l,er0,'->','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.5);
hold on;
semilogy(l,ang0,'-s','Color',[1,0.47,0.1],'MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor',[1.0,0.47,0.1],'LineWidth',1.5);
hold on;
semilogy(l,bnd,'-','Color',[0.4940 0.1840 0.5560],'LineWidth',2);
handle = legend('$\|x_{1}-x_{1}^{(k)}\|/\|R^{-1}\|$',...
    '$\sin\angle(x_{1},x_{1}^{(k)})$','$\kappa(C)\tau/(c_{1}-c_{2})$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
set(gca, 'MinorGridAlpha', 0.02);

figure;
semilogy(l,er0,'->','Color','g','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.5);
hold on;
semilogy(l,er1,'-o','Color','r','MarkerIndices',1:9:k0,...
'MarkerSize',8,'LineWidth',1.5);
hold on;
semilogy(l,er2,'-D','Color','b','MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor','b','LineWidth',1.5);
handle = legend('$\|x_{1}-x_{1}^{(k)}\|/\|R^{-1}\|$', ...
    '$\|x_{1}-\bar{x}_{1}^{(k)}\|/\|R^{-1}\|$',...
    '$\|x_{1}-\tilde{x}_{1}^{(k)}\|/\|R^{-1}\|$','Location','northeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15);
grid on;
set(gca, 'MinorGridAlpha', 0.02);
