function [U, Z, Zb, B, bbeta] = gGKB2_L(A, L, b, k, tol, reorth, type)
    % generalized Golub-Kahan bidiagonalizition (gGKB) of the linear operator:
    %   L: (R(M), <.,.>_M) ---> (R^p, <.,.>_2),  x |--> Lx,
    %   M = A^TA+L^TL maybe positive semidefinite.
    %
    % Inputs:
    %   A, L: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b: right-hand side vector
    %   k: the maximum number of iterations 
    %   tol: stopping tolerance of pcg.m for solving M*sb = s or min||Ms-s||_2
    %       if tol=0, then solve it directly 
    %   reorth: 
    %       0: no reorthogonalization
    %       1: full reorthogonaliation, MGS
    %       2: double reorthogonaliation, MGS
    %   type:
    %       'posi': M is positive definite
    %       'semi': M is positive semidefinite
    %
    % Outputs:
    %   U: mx(k+1) column 2-orthornormal matrix
    %   Z: nx(k+1) matrix, column M-orthonormal 
    %   Zb: nx(k+1) matrix, Zb=MZ
    %   bbeta: 2-norm of b的
    %   B: (k+1)x(k+1) lower bidiagonal matrix
    %
    % Reference: [1]. Haibo Li, Generalizing the Golub-Kahan bidiagonalization for 
    % large-scale GSVD computation.
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 24, March, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 6
        error('Not Enough Inputs')
    end

    if isa(L, 'function_handle') && tol == 0
        error('Tol must not be 0 for a funtional handel L')
    end

    [p, n] = size(L);   [m, ~] = size(A);
    if size(b,1) ~= p || size(A,2) ~= n
        error('The dimensions are not consistent')
    end

    M = A'*A + L'*L;
    if strcmp(type, 'semi') && tol == 0
        Mp = pinv(M);
    end

    % declares the matrix size
    fprintf('[Start gGKB...], max_Iter=%d, reorth=%d\n', [k,reorth]);
    B  = zeros(k+1, k+1);
    U  = zeros(m, k+1);
    Z  = zeros(n, k+1);
    Zb = zeros(n, k+1); 

    % initial step of gGKB
    bbeta = norm(b);
    u = b / bbeta;  
    U(:,1) = u;
    
    rb = L' * u;  
    
    if tol == 0 && strcmp(type, 'posi')
        r = M \ rb;
    elseif tol == 0 && strcmp(type, 'semi')
        r = Mp * rb;
    else
        r = pcg(M, rb, tol, 2*n);
        % r = lsqr(M, rb,tol, 2*n);
    end

    alpha = sqrt(r'*M*r);
    z  = r / alpha;     Z(:,1) = z;
    Zb(:,1) = M * z;
    B(1,1)  = alpha;

    % k step iteration of gGKB
    for j = 1:k
        fprintf('[gGKB iterating...], step=%d--------\n', j);
        % compute u in 2-inner product
        s = L * z - alpha * u;
        if reorth == 1
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
        elseif reorth == 2
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
        end

        beta = norm(s);
        if beta < 1e-14
            fprintf('[Breakdown...], beta=%f, gGKB breakdown at %d--\n', [beta,j]);
            U  = U(:,1:j);
            Z  = Z(:,1:j);
            Zb = Zb(:,1:j);
            B  = B(1:j,1:j);
            break;
        end

        u = s / beta;
        U(:,j+1) = u;
        B(j+1,j) = beta;

        % compute z in M-inner product
        rb = L' * u;
        if tol == 0 && strcmp(type, 'posi')
            r = M \ rb;
        elseif tol == 0 && strcmp(type, 'semi')
            r = Mp * rb;
        else
            r = pcg(M, rb, tol, 2*n);
            % r = lsqr(M, rb,tol, 2*n);
        end

        r = r - beta * z;

        if reorth == 1
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
        elseif reorth == 2
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
        end

        alpha = sqrt(r'*M*r);
        if alpha < 1e-14
            fprintf('[Breakdown...], alpha=%f, gGKB breakdown at %d--\n', [alpha,j]);
            U  = U(:,1:j+1);
            Z  = Z(:,1:j);
            Zb = Zb(:,1:j);
            B  = B(1:j+1,1:j);
            break;
        end

        z = r / alpha;
        Z(:,j+1) = z;
        Zb(:,j+1)  = M * z;
        B(j+1,j+1) = alpha;
    end

end
