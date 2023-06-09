%************************************************************************

%Piecewise Linear Interpolation

%************************************************************************

% Known points
x = [1, 2, 4, 8];
y = [3, 7, 21, 73];

% Query points
query_points = [3, 7];

% Perform piecewise linear interpolation
interpolated_values = interp1(x, y, query_points, 'linear');

% Display results
disp(interpolated_values);

%************************************************************************
%Ans   14    60
%************************************************************************

%Cubic Spline Interpolation

%************************************************************************

x = [0, 1, 2, 3];
y = [1, 2, 33, 244];

% Query points
query_points = 2.5;

% Perform cubic spline interpolation
interpolated_values = spline(x, y, query_points);

% Display results
disp(interpolated_values);

%************************************************************************
%Ans   106.6250
%************************************************************************

%Sterling's  Interpolation

%************************************************************************

% Known points
x = [1, 1.3, 1.6, 1.9, 2.2];
y = [0.76, 0.62, 0.45, 0.28, 0.11];

% Query point
query_point = 1.5;

% Compute divided differences
div_diff = divided_differences(x, y);

% Perform Stirling's interpolation
interpolated_value = stirling_interpolation(x, y, div_diff, query_point);

% Display result
disp(interpolated_value);

% Function to compute divided differences
function div_diff = divided_differences(x, y)
    n = length(x);
    div_diff = zeros(n, n);
    div_diff(:, 1) = y;
    for j = 2:n
        for i = 1:n-j+1
            div_diff(i, j) = (div_diff(i+1, j-1) - div_diff(i, j-1)) / (x(i+j-1) - x(i));
        end
    end
end

% Function to perform Stirling's interpolation
function interpolated_value = stirling_interpolation(x, y, div_diff, query_point)
    n = length(x);
    interpolated_value = y(1);
    term = 1;
    for j = 1:n-1
        term = term * (query_point - x(j));
        interpolated_value = interpolated_value + term * div_diff(1, j+1);
    end
end

%************************************************************************
%Ans   0.5075
%************************************************************************

%Bessel%s Interpolation

%************************************************************************

% Known points
x = [1, 1.5, 2, 2.5, 3.0];
y = [10.24, 12.3452, 15.231, 17.5412, 19.3499];

% Query point
query_point = 1.55;

% Order of Bessel function
v = 1;

% Compute Bessel coefficients using least squares approximation
n = length(x);
A = zeros(n, v+1);
for k = 1:v+1
    A(:, k) = besselj(k-1, x)';
end
coefficients = A \ y';

% Compute central difference coefficients
h = x(2) - x(1);
central_diff_coeffs = [1/12, -2/3, 0, 2/3, -1/12] / h^2;

% Perform Bessel's central difference interpolation
interpolated_value = 0;
for k = 1:v+1
    interpolated_value = interpolated_value + coefficients(k) * besselj(k-1, query_point);
    interpolated_value = interpolated_value + (h^2/2) * central_diff_coeffs(k) * besselj(k-1, query_point);
end

% Display result
disp(interpolated_value);

%************************************************************************
%Ans   14.4109
%************************************************************************
