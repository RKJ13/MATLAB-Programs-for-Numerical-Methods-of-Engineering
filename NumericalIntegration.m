%************************************************************************

%Boole%s formula

%************************************************************************

% Function handle for the integrand
fun = @(x) 1./(1+x^2);

% Integration limits
a = 0; % Lower limit
b = 1; % Upper limit

% Number of subintervals (must be a multiple of 4)
n = 4;

% Step size
h = (b - a) / n;

% Compute the integral using Boole's formulas
integral_approximation = 2 * h / 45 * (7 * sum(fun(a + (0:4:n-4) * h)) + 32 * sum(fun(a + (1:4:n-3) * h)) + 12 * sum(fun(a + (2:4:n-2) * h)) + 32 * sum(fun(a + (3:4:n-1) * h)) + 7 * sum(fun(a + (4:4:n) * h)));

% Display the result
disp(integral_approximation);

%************************************************************************
%Ans   0.7855
%************************************************************************

%Romberg’s integration method with Trapezoidal rule

%************************************************************************

% Function handle for the integrand
fun = @(x) 1./(1+x);

% Integration limits
a = 0; % Lower limit
b = 1; % Upper limit

% Number of iterations
n = 4;

% Initialize the Romberg table
R = zeros(n);

% Compute the trapezoidal rule approximation for the first row
h = b - a;
R(1, 1) = h / 2 * (fun(a) + fun(b));

% Perform Richardson extrapolation
for i = 2:n
    h = h / 2;
    
    % Compute the trapezoidal rule approximation for the current row
    R(i, 1) = 0.5 * R(i-1, 1) + h * sum(fun(a + (1:2^(i-2))*(2*h)));
    
    % Perform Richardson extrapolation for each column
    for j = 2:i
        R(i, j) = R(i, j-1) + (R(i, j-1) - R(i-1, j-1)) / (4^(j-1) - 1);
    end
end

% Approximated integral
integral_approximation = R(n, n);

% Display the result
disp(integral_approximation);


% Approximated integral
integral_approximation = R(n, n);

% Display the result
disp(integral_approximation);

%************************************************************************
%Ans   0.6267

 %Ans   0.6267
%************************************************************************

%Romberg Integration with Simpsons 1/3rd rule

%************************************************************************

% Function handle for the integrand
fun = @(x) 1./(1+x.^2);

% Integration limits
a = 0; % Lower limit
b = 1; % Upper limit

% Number of iterations
n = 1;

% Initialize the Romberg table
R = zeros(n);

% Compute the Simpson's 1/3 rule approximation for the first row
h = (b - a) / 2;
x = a:h:b;
y = fun(x);
R(1, 1) = h / 3 * (y(1) + 4 * sum(y(2:2:end-1)) + 2 * sum(y(3:2:end-2)) + y(end));

% Perform Richardson extrapolation
for i = 2:n
    h = h / 2;
    x = a:h:b;
    y = fun(x);
    
    % Compute the Simpson's 1/3 rule approximation for the current row
    R(i, 1) = 0.5 * R(i-1, 1) + h / 2 * sum(y(1:end-1) + y(2:end));
    
    % Perform Richardson extrapolation for each column
    for j = 2:i
        R(i, j) = R(i, j-1) + (R(i, j-1) - R(i-1, j-1)) / ((2^(2*j-2)) - 1);
    end
end

% Approximated integral
integral_approximation = R(n, n);

% Display the result
disp(integral_approximation);

%************************************************************************
%Ans   0.7833
%************************************************************************

%Double integral using Simpson’s rule 

%************************************************************************

% Integration limits
a = 1;
b = 1.5;
c = 1;
d = 2;

% Step sizes
h = 0.5;
k = 0.25;

% Number of intervals
n = (b - a) / h;
m = (d - c) / k;

% Initialize the integral value
integral_value = 0;

% Compute the double integral using Simpson's rule
for i = 1:n
    for j = 1:m
        x0 = a + (i - 1) * h;
        x1 = a + i * h;
        y0 = c + (j - 1) * k;
        y1 = c + j * k;
        
        integral_value = integral_value + (h * k / 9) * (f(x0, y0) + 4 * f(x0, (y0 + y1) / 2) + f(x0, y1) + ...
            4 * f((x0 + x1) / 2, y0) + 16 * f((x0 + x1) / 2, (y0 + y1) / 2) + 4 * f((x0 + x1) / 2, y1) + ...
            f(x1, y0) + 4 * f(x1, (y0 + y1) / 2) + f(x1, y1));
    end
end

% Display the result
disp(integral_value);

% Function f(x, y) = 1/(x + y)
function value = f(x, y)
    value = 1 / (x + y);
end

%************************************************************************
%Ans   0.7376
%************************************************************************

%Double integral using Trapizoidal rule 

%************************************************************************

% Integration limits
a = 1;
b = 2;
c = 1;
d = 2;

% Step sizes
h = 0.5;
k = 0.5;

% Number of intervals
n = (b - a) / h;
m = (d - c) / k;

% Initialize the integral value
integral_value = 0;

% Compute the double integral using the Trapezoidal rule
for i = 1:n
    for j = 1:m
        x0 = a + (i - 1) * h;
        x1 = a + i * h;
        y0 = c + (j - 1) * k;
        y1 = c + j * k;
        
        integral_value = integral_value + (h * k / 4) * (f(x0, y0) + f(x0, y1) + f(x1, y0) + f(x1, y1));
    end
end

% Display the result
disp(integral_value);

% Function f(x, y) - Define your integrand function here
function value = f(x, y)
    % Example function: f(x, y) = x + y
    value = x + y;
end

%************************************************************************
%Ans   3
%************************************************************************