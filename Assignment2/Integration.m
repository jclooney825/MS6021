%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of subintervals
N = 10000;

% Interval values: [a,b] 
a = 0;
b = 1; 

% Value of exponent value
alpha = 2.2;

% Function handle 
funct = @(x) x.^(alpha).*(1 + 3*x.^2); 

% Integral value 
I_alpha = trapozodial(funct, a, b, N)

% Resulting error
err = errors(funct, a, b, N, alpha)

% Trapozodial Integration function 
function result = trapozodial(funct, a, b, N)
    format long g 

    % Step size 
    delta_x = (b - a)/N; 

    % x values from bounds a to b
    x = linspace(a,b, N+1);
    
    % Generated values 
    vals = 2*funct(x);

    % Multiply bound values by 1/2 
    vals(1) = (1/2)*(vals(1)); 
    vals(end) = (1/2)*(vals(end)); 

    % Sums the resulting values and returns them 
    result = (delta_x/2)*sum(vals);

end 

% Exact solution of our integral 
function I = exact_solution(alpha)
    I = ((1)/(alpha + 1)) + ((3)/(alpha + 3));  
end 

% Function that finds the error
function err = errors(funct, a, b, N, alpha)
    % Analytical and numerical solutions 
    I_exact = exact_solution(alpha);
    I_numerical = trapozodial(funct, a, b, N);
   
    % Error = max(|computed solution - exact solution|)
    err = abs(I_exact - I_numerical);
end 
