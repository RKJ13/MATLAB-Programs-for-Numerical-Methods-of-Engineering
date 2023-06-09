%**************************************************************************

%Richardson’s extrapolation

%**************************************************************************

% Function handle
fun = @(x) 5 * x .* exp(-2 * x);

% Query point
query_point = 0.35;

% Step sizes
h = [ 0.125, 0.25];

% Compute central difference approximations
approximations = (fun(query_point + h) - fun(query_point - h)) ./ (2 * h);

% Perform Richardson's extrapolation if enough elements are available
if numel(approximations) >= 2
    extrapolated_approximation = (4 * approximations(2) - approximations(1)) ./ 3;
    
    % Display result
    disp(extrapolated_approximation);
else
    disp("Not enough elements in the approximations array to perform Richardson's extrapolation.");
end

%**************************************************************************
%Ans   1.0497
%************************************************************************
