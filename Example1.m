% show the error \tilde{g}_i of recurrence of v_til
% and loss of orthogonality of v_til;
% JBD do not use any reorthogonalization
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

% N = 10000;
% c = N:-1:1;  c = c/(2*N);
% s= sqrt(1 - c.*c);
% C = diag(c);  S = diag(s);  
% kapp = 5;
% RR = diag(linspace(1,kapp,N));
% invR = diag(1./linspace(1,kapp,N));
% AA = C*RR;  LL = S*RR;
% A = sparse(AA);  L1 = sparse(LL);
% X = invR; 


%%---------JBD process------------------------------------------------
L=0.1*L1;
[m,n]=size(A);  p=size(L,1);
C = [A;L];  
[Q,R]=qr1(C);  kap = cond(R);
% QQ = [C; S]; Q = sparse(QQ);  kap = kapp;
% Qa = sparse(C);  Ql = sparse(S);

k0=60;
tol1 = 1e-12;
tol2 = 1e-9;
reorth=0;
rng(2022);  % random seed
t = randn(n, 1);

[B, Bh, U, Uh, V_til]=JBD(A,L,t,k0+1,tol1,reorth);
[B2, Bh2, U2, Uh2, V2_til]=JBD(A,L,t,k0+1,tol2,reorth);
P=eye(k0);
for i=1:k0 
    P(i+1,i+1)=-1*P(i,i);
end

s1 = zeros(k0, 1);  % error of recurrences of v_til
s11 = zeros(k0, 1);  % error of recurrences of v1_til
s2 = zeros(k0, 1);  % orthogonality level of V_til
s22 = zeros(k0, 1);  % orthogonality level of V2_til
Q_p = Q*Q';
for i=1:k0
    err = Q_p*[U(:,i);zeros(p,1)] - B(i,i)*V_til(:,i) - B(i,i+1)*V_til(:,i+1);
    s1(i) = norm(err);
    err2 = Q_p*[U2(:,i);zeros(p,1)] - B2(i,i)*V2_til(:,i) - B2(i,i+1)*V2_til(:,i+1);
    s11(i) = norm(err2);
    H = eye(i) - V_til(:,1:i)'* V_til(:,1:i);
    s2(i) = norm(H);
    H2 = eye(i) - V2_til(:,1:i)'* V2_til(:,1:i);
    s22(i) = norm(H2);
end

%%-----------plot--------------------------------------
lw = 2; l = 1:k0;
figure; 
semilogy(l,s1,'-.b','LineWidth',2);
hold on;
semilogy(l,3*kap*tol1*ones(k0,1),'-r','LineWidth',2)
semilogy(l,s11,'-.','Color',[0.3010 0.7450 0.9330],'LineWidth',2);
hold on; %[0.4940 0.1840 0.5560]
semilogy(l,3*kap*tol2*ones(k0,1),'-','Color',[0.9500 0.6 0],'LineWidth',2)
handle=legend('$\|\tilde{g}_{i}\|$, $\tau_{1}$','$3\kappa(C)\tau_{1}$',...
    '$\|\tilde{g}_{i}\|$, $\tau_{2}$','$3\kappa(C)\tau_{2}$','Location','southeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Error','Fontsize',15); 
grid on;
set(gca, 'MinorGridAlpha', 0.02);

figure; 
semilogy(l,s2,'-s','Color',[0.6350 0.0780 0.1840],'MarkerIndices',1:9:k0,...
'MarkerSize',8,'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',1.5);
hold on; 
semilogy(l,s22,'-o','Color',[0 0 0.8],'MarkerIndices',1:9:k0,...
    'MarkerSize',6,'MarkerFaceColor',[0 0 0.8],'LineWidth',1.5);
handle=legend('$\|I_{k}-\widetilde{V}_{k}^{T}\widetilde{V}_{k}\|$, $\tau_{1}$',...
    '$\|I_{k}-\widetilde{V}_{k}^{T}\widetilde{V}_{k}\|$, $\tau_{2}$','Location','southeast');
set(handle,'Fontsize',14,'interpreter','latex');
xlabel('Iteration','Fontsize',15);
ylabel('Orthogonality level','Fontsize',15);
grid on;
set(gca, 'MinorGridAlpha', 0.02);

