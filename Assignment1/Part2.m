%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical Modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of selected N values 
N_list = [8,16,32,64,128,256,512];

plot_all(N_list)
plot_errs(N_list)

% Main numerical solution function
function [x,U] = numerical_result(N)

    % Step size 
    h = 1/N;

    % Boundary conditions 
    u1 = 0;
    u_end = 0; 
    
    % iteration values 
    x_i = linspace(h,1-h,N-1);
    
    % a(x) = exp(x^2 - x) 
    diag = 2 + h^2*exp(x_i.^2 - x_i)';

    % Constructing a sparse matrix 
    e = ones(N,1);
    A = spdiags([-e e -e],-1:1,N-1,N-1);
    d = 0;
    A = spdiags(diag,d,A); 
  
    % F column vector values
    % f(x) = 7/(2 + x^2 + x) 
    F = h^2*(7./(2 + x_i.^2 + x_i))';
    
    % Solution of matrix equation 
    U = [inv(A)*F]';
    
    % x values
    x = linspace(0,1,N+1);
    
    % Add in boundary values
    U = [u1,U,u_end]';
 
end 

% Double mesh error calculation 
function double_mesh_err = double_mesh_error(N)
    format long g
    
    % Numerical solutions for N and 2N values
    [x,U] = numerical_result(N);
    [x,U_2N] = numerical_result(2*N);
    
    % Uses every second value of 2N solution 
    U_2N = U_2N(1:2:end);
    
    % Max double mesh error
    double_mesh_err = max(abs(U - U_2N)); 
    
end 

% Plots errors as N increases 
function plot_errs(N_list)
    
    % Error vector 
    errs = arrayfun(@double_mesh_error,N_list)'; 
    f1 = figure;
    semilogy(N_list,errs)
    
    % Plot settings 
    xlim([N_list(1), N_list(end)]);
    xlabel('N', 'FontSize',15, 'Interpreter','latex');
    ylabel('log(DME)', 'FontSize',15, 'Interpreter','latex');
    title('Double Mesh Error increasing N','FontSize', 18,'Interpreter','latex');
    grid on;
    print('-depsc2','-painters','figures/errors_part2.eps'); 
end 

% Plots all results in N_list
function plot_all(N_list)

    hold on
    legend('FontSize',12, 'Interpreter','latex')
    % Goes through list and plots for different N
    for i = N_list 
        [x, u] = numerical_result(i);
        plot(x,u, '-','DisplayName', strcat('N=', num2str(i)));
    end 
    
    % Plot settings 
    xlabel('$x$', 'FontSize',15, 'Interpreter','latex') 
    ylabel('$u(x)$','FontSize',15,'Interpreter','latex')
    title(strcat('Numerical Solution of BVP with increasing values of N', ''),'FontSize',18, 'Interpreter','latex')
    grid on 
    print('-depsc2','-painters','figures/part2_sol.eps'); 
    hold off 
    
end 