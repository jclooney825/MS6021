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

% Calling the function 
[t,y] =  Runge_Kutta(f, tspan, y0, N); 

% Plot results 
plot(t,y)

% Plot settings
grid on 
hold on 
xlabel('$t$','Interpreter','latex','FontSize', 15)
ylabel('$y$','Interpreter','latex','FontSize', 15)
title('Problem 1: 4th Order Runge Kutta Method with $N = 600$','Interpreter','latex', 'FontSize', 18)

function [t,y] = Runge_Kutta(f, tspan, y0, N)
    
  % Step size 
    h = (tspan(2) - tspan(1))/N;
    
    % Set first elements of t and y 
    t = linspace(tspan(1),tspan(2), N+1); 
    y(:,1) = y0; 
   
    % 4th Order Runge Kutta Method 
    for i = 1:N
        k1 = f(t(i),y(:,i));
        k2 = f(t(i) + h/2, y(:,i) + h*(k1/2)); 
        k3 = f(t(i) + h/2, y(:,i) + h*(k2/2)); 
        k4 = f(t(i) + h, y(:,i) + h*k3); 
    
        y(:,i+1) = y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    t = t';
    y = y';
end 



