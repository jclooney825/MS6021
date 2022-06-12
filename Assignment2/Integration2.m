%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of subintervals
N = 10000;

% g(x) function 
g = @(x) cos(5 + x.^3); 

% Value of alpha 
alpha = -0.9; 

% Integral value
I_g = tailored(g, N, alpha)

function result = tailored(g, N, alpha)
    format long 
    % Interval values: [a,b] 
    a = 0; 
    b = 1;

    % x values
    x = linspace(a,b,N);
    
    % x_(i+1) and x_(i) values 
    x_i = x(1:N-1);
    x_iplus1 = x(2:N);

    % Coefficient of interpolated line: g^I(x) = Ax + B
    A = (g(x_iplus1) - g(x_i))./(x_iplus1 - x_i);
    B =  g(x_i) - A.*x_i; 
    
    % Integrals over region x_(i) and x_(i+1) 
    first_integral = (A/(alpha + 2)).*((x_iplus1).^(alpha + 2) - (x_i).^(alpha + 2));
    second_integral = (B/(alpha + 1)).*((x_iplus1).^(alpha + 1) - (x_i).^(alpha + 1));

    % Sum over all intervals 
    result = sum(first_integral + second_integral);

end 

