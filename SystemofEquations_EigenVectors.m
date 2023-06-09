%************************************************************************

% Tridiagonal systems using Thomas algorithm

%************************************************************************

a = [0, -0.4, -0.4];  % Lower diagonal
b = [0.8, 0.8, 0.8];  % Main diagonal
c = [-0.4, -0.4, 0];  % Upper diagonal
d = [41, 25, 105];  % Right-hand side

solution = thomas_algorithm(a, b, c, d);
disp(solution);

function x = thomas_algorithm(a, b, c, d)
    n = length(d);
    c_dash = zeros(n, 1);
    d_dash = zeros(n, 1);
    x = zeros(n, 1);

    % Step 1: Forward elimination
    c_dash(1) = c(1) / b(1);
    d_dash(1) = d(1) / b(1);

    for i = 2:n
        temp = b(i) - (a(i) * c_dash(i - 1));
        c_dash(i) = c(i) / temp;
        d_dash(i) = (d(i) - (a(i) * d_dash(i - 1))) / temp;
    end

    % Step 2: Backward substitution
    x(n) = d_dash(n);

    for i = n-1:-1:1
        x(i) = d_dash(i) - (c_dash(i) * x(i + 1));
    end
end


%************************************************************************
%Ans  173.7500
%Ans  245.0000
%Ans  253.7500
%************************************************************************

%non-linear equations using Newton’s method

%************************************************************************

% Define the system of equations
f = @(x) [
    x(1)^2 + x(2)*x(1) + x(2)^2 - 7;
    x(1)^3 + x(2)^3 - 9
];

% Define the Jacobian matrix
J = @(x) [
    2*x(1)+x(2), 2*x(2)+x(1);
    3*x(1)^2, 3*x(2)^2
];

% Initial guess
x0 = [1.5; 0.5];

% Call the Newton's method function
x = newtons_method(f, J, x0, 100, 1e-6);

function x = newtons_method(f, J, x0, max_iterations, tolerance)
    % f: Function handle for the system of equations
    % J: Function handle for the Jacobian matrix
    % x0: Initial guess
    % max_iterations: Maximum number of iterations
    % tolerance: Tolerance for convergence
    
    x = x0;
    iteration = 0;
    error = Inf;
    
    while error > tolerance && iteration < max_iterations
        delta_x = -J(x) \ f(x);
        x = x + delta_x;
        
        error = norm(delta_x);
        iteration = iteration + 1;
    end
    
    if iteration >= max_iterations
        disp("Maximum iterations reached without convergence.");
    else
        disp("Converged to a solution:");
    end
    
    disp(x);
end

%************************************************************************
%Ans  2
%Ans  1
%************************************************************************

%Fixed Point Iteration Method

%************************************************************************

% Define the fixed-point iteration function
g = @(x) 1 ./ sqrt(x+1);

% Initial guess
x0 = 0;

% Call the Fixed Point Iteration Method function
x = fixed_point_iteration(g, x0, 100, 1e-6);
function x = fixed_point_iteration(g, x0, max_iterations, tolerance)
    % g: Function handle for the fixed-point iteration function
    % x0: Initial guess
    % max_iterations: Maximum number of iterations
    % tolerance: Tolerance for convergence
    
    x = x0;
    iteration = 0;
    error = Inf;
    
    while error > tolerance && iteration < max_iterations
        x_new = g(x);
        
        error = abs(x_new - x);
        x = x_new;
        
        iteration = iteration + 1;
    end
    
    if iteration >= max_iterations
        disp("Maximum iterations reached without convergence.");
    else
        disp("Converged to a solution:");
    end
    
    disp(x);
end
%************************************************************************
%Ans  0.7549
%************************************************************************

%Rayleigh Power method

%************************************************************************

% Define the matrix
A = [5, 4, 2;4, 5, 2;2, 2, 2;];

% Initial guess for the dominant eigenvector
x0 = [1; 0; 0;];

% Call the Rayleigh Power Method function
[lambda, v] = rayleigh_power_method(A, x0, 100, 1e-6);

function [lambda, v] = rayleigh_power_method(A, x0, max_iterations, tolerance)
    % A: Square matrix
    % x0: Initial guess for the dominant eigenvector
    % max_iterations: Maximum number of iterations
    % tolerance: Tolerance for convergence
    
    x = x0;
    iteration = 0;
    error = Inf;
    
    while error > tolerance && iteration < max_iterations
        x_new = A * x;
        x_new = x_new / norm(x_new);
        
        lambda = x_new' * A * x_new;
        
        error = norm(x_new - x);
        x = x_new;
        
        iteration = iteration + 1;
    end
    
    v = x;
    
    if iteration >= max_iterations
        disp("Maximum iterations reached without convergence.");
    else
        disp("Converged to an eigenvalue and eigenvector pair:");
    end
    
    disp("Eigenvalue:");
    disp(lambda);
    
    disp("Eigenvector:");
    disp(v);
end

%************************************************************************
%{
Ans :
    Eigenvalue:
    10.0000

    Eigenvector:
    0.6667
    0.6667
    0.3333
%}
%************************************************************************

%Jacobi%s Method

%************************************************************************

% Define the symmetric matrix A
A = [2, 1, 1; 1, 2, 1; 1, 1, 2];

% Call the Jacobi eigenvalue algorithm function
[eigenvalues, eigenvectors] = jacobi_eigenvalue_algorithm(A, 1000, 1e-6);

function [eigenvalues, eigenvectors] = jacobi_eigenvalue_algorithm(A, max_iterations, tolerance)
    % A: Symmetric matrix
    % max_iterations: Maximum number of iterations
    % tolerance: Tolerance for convergence
    
    n = size(A, 1);
    D = A; % Initialize the diagonal matrix
    V = eye(n); % Initialize the eigenvector matrix
    
    iteration = 0;
    error = Inf;
    
    while error > tolerance && iteration < max_iterations
        max_off_diag = max(max(triu(abs(D), 1))); % Find maximum off-diagonal element
        
        if max_off_diag <= tolerance
            break; % Convergence condition
        end
        
        [p, q] = find(abs(D) == max_off_diag, 1); % Find indices of maximum off-diagonal element
        
        theta = atan2(2 * D(p, q), D(q, q) - D(p, p)) / 2; % Compute rotation angle
        
        G = eye(n);
        G(p, p) = cos(theta);
        G(q, q) = cos(theta);
        G(p, q) = sin(theta);
        G(q, p) = -sin(theta); % Compute Givens rotation matrix
        
        D = G' * D * G; % Update diagonal matrix
        V = V * G; % Update eigenvector matrix
        
        error = max_off_diag;
        iteration = iteration + 1;
    end
    
    eigenvalues = diag(D); % Extract eigenvalues from diagonal matrix
    eigenvectors = V; % Eigenvectors are stored in the eigenvector matrix
    
    if iteration >= max_iterations
        disp("Maximum iterations reached without convergence.");
    else
        disp("Converged to eigenvalues and eigenvectors:");
    end
    
    disp("Eigenvalues:");
    disp(eigenvalues);
    
    disp("Eigenvectors:");
    disp(eigenvectors);
end

%************************************************************************
%{
Ans :
  Eigenvalues:
     4
     1
     1

Eigenvectors:
    0.5774   -0.7071   -0.4082
    0.5774    0.7071   -0.4082
    0.5774         0    0.8165

%}
%************************************************************************

%Given’s method

%************************************************************************

% Define the input matrix
A = [2, 1, -2; 1, 2, -2; -2, -2, 3];

% Call the givens_tridiagonal function
T = givens_tridiagonal(A);

% Display the tridiagonal form
disp("Tridiagonal Form:");
disp(T);
function T = givens_tridiagonal(A)
    % A: Input matrix
    
    [m, n] = size(A);
    
    if m ~= n
        error("Input matrix must be square.");
    end
    
    T = A;
    
    for k = 1:n-2
        for i = k+2:n
            if T(i, k) ~= 0
                r = hypot(T(k+1, k), T(i, k));
                c = T(k+1, k) / r;
                s = -T(i, k) / r;
                
                T([k+1, i], :) = [c, -s; s, c] * T([k+1, i], :);
                T(:, [k+1, i]) = T(:, [k+1, i]) * [c, -s; s, c]';
            end
        end
    end
end

%************************************************************************
%{
Ans :
    Tridiagonal Form:
    2.0000    2.2361         0
    2.2361    4.4000    0.8000
         0    0.8000    0.6000
  

%}
%************************************************************************
