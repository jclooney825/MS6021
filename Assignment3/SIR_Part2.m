%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for Part 2 
data = [
40,1,0,12,2,518,2.9;
41,0,2,8,0,529,1.89;
42,1,0,16,1,612,2.94;
43,0,0,15,2,656,2.59;
44,1,0,24,5,746,4.02;
45,0,2,10,3,578,2.6;
46,1,1,15,0,666,2.55;
47,0,0,18,5,617,3.73;
48,0,3,17,2,636,3.46;
49,0,0,17,6,806,2.85;
50,1,1,1,4,686,1.02;
51,8,3,0,11,904,2.43;
52,3,0,45,17,1037,6.27;
53,8,5,41,18,1230,5.85;
54,4,13,66,18,1164,8.68;
55,17,32,97,18,1291,12.7;
56,29,46,183,43,1545,19.48;
57,33,133,598,83,3194,26.52;
58,64,293,726,112,3651,32.73;
59,45,271,780,157,3752,33.4;
60,56,354,735,254,4159,33.64;
61,33,238,537,234,3401,30.64;
62,16,194,362,202,3001,25.79;
63,6,75,166,151,1923,20.7;
64,21,119,115,114,1903,19.39;
65,2,32,48,126,1157,17.98;
66,1,1,44,72,925,12.76;
67,0,4,25,40,842,8.19;
68,0,3,31,28,681,9.1;
69,0,0,18,15,569,5.8;
70,1,1,19,17,526,7.22;
71,0,0,2,3,371,1.35;
72,0,0,1,4,296,1.69]; 

plot_results(data)

function [R_best, M] = find_best(t0, data)
    
    % Population size 
    population = 5.2*10^(7);
    
    % Average recovery period 
    inv_gamma = 3; % days; 

    R_number = [1.21:0.001:1.24]';
    
    Beta = (1/inv_gamma)*R_number;

    % Time for disease spread 
    tspan = [t0, 7*(46:72)];

    % Initial conditions 
    I0 = 1/population;
    y0 = [1 - I0, I0]; 

    % Length of R_number vector 
    k = max(size(R_number));

    % A row of repeated intial conditions for our ODE45 function
    y0 = repelem(y0,k);
    
    % Function handle 
    f = @(t,y) [-Beta.*y(1:k).*y(k+1:2*k); Beta.*y(1:k).*y(k+1:2*k) - (1/inv_gamma)*y(k+1:2*k)];    

    % Ode45 solverolutions 
    options = odeset('RelTol',1e-13,'AbsTol',1e-20); 
    [t,y] = ode45(f, tspan, y0, options); 

    % Selecting week number and weekly incidence from our data matrix 
    week_number = data(:,1);
    weekly_incidence =  data(:,5);
    
    % Removing unwanted data 
    week_number = week_number(8:end);
    weekly_incidence =  weekly_incidence(8:end);
    
    % Incdidence modeled with data discarded at t0. 
    S = y(:,1:k);
    incidence_modelled = -diff(S(2:end,1:k));
    
    % Normalization factor 
    norm = sum(weekly_incidence)./sum(incidence_modelled);
    
    % Normalized incidence modelled
    incidence_modelled = norm.*incidence_modelled;
    
    % Return minimum value of least squares 
    sum_of_squares = sum((weekly_incidence - incidence_modelled(:,1:k)).^2);
    [M,I] = min(sum_of_squares);
   
    % Return best R value 
    R_best = R_number(I);
end 

% Finds the best values of R and t0 
function [t_best, R_best] = find_best_t0(data) 
    val = inf; 
    for t0 = 200:1:250
        [R_val, least_squares] = find_best(t0, data);
        if least_squares < val
            t_best = t0;
            R_best = R_val;
        end 
        val = least_squares;
    end
end 

function plot_results(data)

    [t_best, R_best] = find_best_t0(data);
    % Population size 
    population = 5.2*10^(7);
    
    % R number 
    R_number = R_best;
    
    % Average recovery period 
    inv_gamma = 3; % days; 
    
    Beta = (1/inv_gamma)*R_number;
    
    % Inital start date 
    t0 = t_best;
    
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

    % Selecting week number and weekly incidence from our data matrix 
    week_number = data(:,1);
    weekly_incidence =  data(:,5);
    
    % Removing unwanted data 
    week_number = week_number(8:end);
    weekly_incidence =  weekly_incidence(8:end);
    
    f2 = figure; 
    
    % Incdidence modeled with data discarded at t0. 
    incidence_modelled = -diff(S(2:end));
    
    % Normalization factor 
    norm = sum(weekly_incidence)/sum(incidence_modelled);
    
    % Normalized incidence modelled
    incidence_modelled = norm*incidence_modelled;
    
    % Areas under data and model plots 
    area = 7*sum(weekly_incidence);
    area_modeled = 7*sum(incidence_modelled);
    
    sum_of_squares = sum((weekly_incidence - incidence_modelled).^2);
    
    % Plotting results 
    hold on 
    grid on 
    plot(7*week_number, weekly_incidence, '-o', 'Color', 'black')
    plot(7*week_number,incidence_modelled, 'Color', 'red')
    
    xlabel('Time (days)','Interpreter','latex', 'FontSize', 15)
    ylabel('Cases','Interpreter','latex', 'FontSize', 15)
    title('SIR Model for B influenza with Fitted R and $t_0$','Interpreter','latex', 'FontSize', 18)
    legend('Data', 'Model','Interpreter','latex','FontSize', 12);
     
end 
