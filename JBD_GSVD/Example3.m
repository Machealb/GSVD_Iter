% Show the ghost generalized singular values phenomenon by JBD-GSVD algorithm
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences, Dec 05, 2022.

clear;
directory = pwd;
path(directory, path)
path([directory, '/Matrices'], path)

%% {A , L}
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
m = N; n = N; p = N;

%%
L=1*L1;
[Q,R]=qr1([A;L]);
Qa=Q(1:m,:);  Ql=Q(m+1:m+p,:);

k0 = 80;
H1 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at top)
H2 = zeros(k0, k0);   % every column store the approxiamte at the k-th step(store at bottom)

%%  k0-step JBD algorithm
tol = 1e-10;
reorth = 1;
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

for k = 1:k0
    fprintf('Running sub_SVD algorithm: the %d-th step ===================\n', k);
%     B_k = B(1:k,1:k);
%     c_k = svd(B_k);  
%     H1(1:k, k) = c_k;  % approximate singular value of Qa (largest)
%     H2(k0-k+1:k0, k) = c_k;  % approximate singular value of Qa (smallest)

    Bh_k = Bh(1:k,1:k);
    s_k = svd(full(Bh_k));
    %H1(1:k, k) = s_k;  % approximate singular value of Ql (largest)
    H2(k0-k+1:k0, k) = s_k;  % approximate singular value of Ql (smallest)
end

%% plotting
for i = 1: 3
    x_i = i : 1 : k0;
%     plot(x_i, H1(i, i:k0), 'x-');
%     hold on;
    plot(x_i, H2(k0-i+1, i:k0), 'x-');
    hold on;
end
%scatter(k0*ones(N, 1), c, 'ko');
scatter(k0*ones(N, 1), s, 'ko');
xlabel('Iteration','Fontsize',15);
ylabel('Ritz value','Fontsize',15);





