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
[t,y] =  Predictor_Corrector(f, tspan, y0, N);

% Plot results 
plot(t,y(:,1))

% Plot settings 
grid on 
xlabel('$t$','Interpreter','latex','FontSize', 15)
ylabel('$y$','Interpreter','latex','FontSize', 15)
title('Problem 1: Predictor/Corrector Method with $N = 600$','Interpreter','latex', 'FontSize', 18)


function [t,y] = Predictor_Corrector(f, tspan, y0, N)
    
    % Step size 
    h = (tspan(2) - tspan(1))/N;
   
    % Set first elements of t and y 
    t = linspace(tspan(1),tspan(2), N+1); 
    y(:,1) = y0; 

    % Predictor-Corrector method 
    for i = 1:N
        y_star = y(:,i) + h*f(t(i),y(:,i));
        y(:,i+1) = y(:,i) + 0.5*h*( f(t(i),y(:,i)) + f(t(i+1),y_star) );
    end 
    t = t';
    y = y';
end 
