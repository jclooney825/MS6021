%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Population size 
population = 5.2*10^(7);

% R number 
R_number = 1.3;

% Average recovery period 
inv_gamma = 3; % days; 
Beta = (1/inv_gamma)*R_number;

% Inital start date 
t0 = 250; 

% Time for disease spread 
tspan = [t0, 7*(46:72)];

% Initial conditions 
I0 = 1/population;
y0 = [1 - I0, I0]; 

% SIR model functions 
f = @(t,y) [-Beta*y(1)*y(2); Beta*y(1)*y(2) - (1/inv_gamma)*y(2)];

% Matlab's built in ODE solver 
options = odeset('RelTol',1e-13,'AbsTol',1e-20); 
[t,y] = ode45(f, tspan, y0, options);

% Output values for S, I, and R
S = y(:,1);
I = y(:,2);
R = 1 - S - I; 

% Plot results 
f1 = figure;
hold on
plot(t, S, 'Color','blue');
plot(t, I, 'Color','red');
plot(t, R, 'Color','green');
grid on;

% Plot settings
xlim([t0,504])
xlabel('Time (days)', 'FontSize', 15,'Interpreter', 'latex')
title('SIR Model for B influenza in the Midwest (2007-2008)', 'FontSize', 18,'Interpreter', 'latex'); 
legend('Susceptible', 'Infected', 'Recovered','FontSize',12,'Location','northwest','Interpreter', 'latex');
print('-depsc2','-painters','figures/sir_model.eps');