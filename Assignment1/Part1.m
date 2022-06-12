%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical Modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of selected N values 
N_list = [8,16,32,64,128,256,512];

% Calling the functions 
plot_all(N_list)
plot_errs(N_list)

% Main numerical solution function
function [x,U] = numerical_result(N)

    % Step size 
    h = 1/N;

    % Boundary conditions 
    u1 = 0;
    u_end = 0; 
    
    % Iteration values (x_2 to x_N)
    x_i = linspace(h,1-h,N-1);
    
    % Constructing a sparse matrix 
    e = ones(N,1);
    A = spdiags([-e 2*e -e],-1:1,N-1,N-1);
    
    % F column vector values 
    F = [h^(2)*cos(2*x_i)]';
    
    % Solution of matrix equation 
    U = [inv(A)*F]';
    
    % x values 
    x = linspace(0,1,N+1);
    
    % Add in boundary values
    U = [u1,U,u_end]';
    
end 

% Analytical solution 
function [x,y] = exact_solution(N)
    x = linspace(0,1,N+1);
    y = (1/4)*(cos(2*x)+ (1-cos(2))*x - 1);  
end 

% Max error value
function max_err = errors(N)
    
    % Analytical and numerical solutions 
    [x,y] = exact_solution(N);
    [x_n,U] = numerical_result(N);
   
    % Error = max(|computed solution - exact solution|)
    max_err = max(abs(U - y'));
    
end 

% Errors vector as N increases 
function plot_errs(N_list)
    
    % Error plot 
    errs = arrayfun(@errors,N_list)'; 
    f1 = figure;
    semilogy(N_list,errs)
    
    % Error plot settings 
    xlim([N_list(1), N_list(end)])
    xlabel('N', 'FontSize',15, 'Interpreter','latex');
    ylabel('log(Error)', 'FontSize',15, 'Interpreter','latex');
    title('Max Nodal Error with increasing N','FontSize',18,'Interpreter','latex');
    grid on;
end 

% Plots all results in N_list
function plot_all(N_list)

    % Plots the analytical solution first 
    [x,y]= exact_solution(1000);
    hold on 
    plot(x,y, 'black')
    legend('Analytical','FontSize',12, 'Interpreter', 'latex')
    
    % Goes through N_list and plots for different N
    for i = N_list 
        [x, u] = numerical_result(i);
        plot(x,u, '--','DisplayName', strcat('N=', num2str(i)))
    end 
    
    % Plot settings 
    xlabel('x', 'FontSize',15, 'Interpreter','latex') 
    ylabel('u(x)','FontSize',15, 'Interpreter','latex')
    title('Numerical Solution of BVP with increasing values of N','FontSize',18,'Interpreter','latex')
    grid on; 
    hold off 
    
end 


