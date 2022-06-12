%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of subintervals 
N = 600; 

% Initial conditions 
y0 = [1]; 

% Function handle f(t,y) 
% Problem 1 
f = @(t,y) [-y(1)-5*exp(-t)*sin(5*t)]; 

% Time interval
tspan = [0,3]; 

% Calling the Eulers Method function 
[t,y] = Eulers_Method(f, tspan, y0, N);

% Plot solution 
plot(t,y)

% Plot settings
hold on
grid on
xlabel('t','Interpreter','latex','FontSize', 15)
ylabel('$y$','Interpreter','latex','FontSize', 15, 'Rotation',0)
title('Problem 1: Explicit Euler Method with $N = 600$','Interpreter','latex', 'FontSize', 18)

% Eulers Method 
function [t,y] = Eulers_Method(f, tspan, y0, N)

    % Step size 
    h = (tspan(2) - tspan(1))/N;
    
    % Set first elements of t and y 
    t = linspace(tspan(1),tspan(2), N+1); 
    y(:,1) = y0; 

    % For loop that generates t and y vectors 
    for i = 1:N
        y(:,i+1) = y(:,i) + h*f(t(i),y(:,i));
    end 
    % Transpose vectors before returning
    t = t';
    y = y';
end 

