%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Mathematical modeling MSc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size of matrix 
N = 80;

% Plot of eigenvalues for a given N
plot_all_eigs(N)

% Plot the 30 lowest eigenvalues for increadsing N
num_lowest = 30; 
N_list = [31, 60, 300, 3000, 30000];

% Plot the 30 lowest eigs for all numbers in N_list
plot_smallest_30_eigs(N_list, num_lowest);

% A method that returns the Eigenvalues of an N-1 linear operator matrix L 
function eigenvalues = find_eigs(N, num_lowest) 
    % Step size 
    h = 1/N; 
    
    % Constructing a sparse matrix 
    e = ones(N,1);
    L = spdiags([-e/h^2 2*e/h^2 -e/h^2],-1:1,N-1,N-1); 

    % Return the sorted eigenvalues of the discretized linear operator L 

    % Return all eigenvalues 
    if num_lowest == 'all'
        eigenvalues = sort(eig(L));
    % Return smallest num_lowest number of eigs 
    elseif isa(num_lowest, 'float')
        eigenvalues = sort(eigs(L,num_lowest, 'sm')); 
    % Return NaN if any other input is entered
    else
        disp('Input invalid. Please try again.')
        eigenvalues = NaN; 
    end 
end 

% Plots the eigenvalues
function plot_all_eigs(N)
    f1 = figure;
    hold on 
    % Calcualtes Eigenvalues 
    eigenvalues = find_eigs(N, 'all');

    % N vector 
    N_vec = linspace(1,N, N-1);

    % Analytical eigenvalues
    y = N_vec.^2*(pi)^2;

    % Plot both
    plot(N_vec, y, 'o', 'DisplayName','Analytic'); 
    plot(N_vec, eigenvalues, 'o', 'DisplayName', strcat('N = ', num2str(N)));

    % Plot settings 
    xlabel('$n$','FontSize',15, 'Interpreter','latex') 
    ylabel('$\lambda_n$','FontSize',15, 'Interpreter','latex','Rotation',0)
    legend('FontSize',15,'Location','northwest','Interpreter','latex');
    title(strcat('Plot of Eigenvalues with $N = $',string(N)),'FontSize',18,'Interpreter','latex')
    grid on  

end 

% Plots the eigenvalues vs N-1 
function plot_smallest_30_eigs(N_list, num_lowest)
    
    % Create figure
    f2 = figure;
    hold on; 
    legend('FontSize',12,'Location','northwest', 'Interpreter','latex');
   
    % Vector on integers 
    n_vec = linspace(1,num_lowest, num_lowest);  
    
    % Analytical eigenvalues
    analytical_eigs = n_vec.^2*(pi)^2; 
    
    % Plot analytical solution and label on legend 
    plot(n_vec, analytical_eigs, 'o', 'Color','black');
    legend('Analytic','FontSize',12);

    % Calcualtes lowest Eigenvalues in N_list 
    for i = N_list
        % Loops through and plots 30 lowest of
        % N eigenvalues for increasing values 
        eigenvalues = find_eigs(i, num_lowest);
        plot(n_vec, eigenvalues, 'o','DisplayName', strcat('N=', num2str(i)));
    end

    % Plot settings 
    xlabel('$n$','FontSize',15, 'Interpreter','latex');
    ylabel('$\lambda_n$','FontSize',15, 'Interpreter','latex', 'Rotation',0);
    title('$30$ Smallest Eigenvalues with increasing $N$','FontSize',18, 'Interpreter','latex');
    grid on

end 

