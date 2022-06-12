%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% James Clooney 
% Questions 6 - 8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MS6021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size of matrix 
N = 80;

% Number of smallest eigenvalues used
num_lowest = 30; 

% List of N values
N_list = [31, 60, 300, 3000, 30000];

% Plots first thirty eigenvalues for increasing N
plot_smallest_30_eigs(N_list, num_lowest)

% A method that returns the Eigenvalues of an N-1 linear operator matrix L 
function eigenvalues = find_eigs(N, num_lowest)   
    % Step size 
    h = 1/N; 
    
    % Iteration values
    x_i = linspace(h,1-h,N-1);

    % Diagonal values of our matrix 
    diag = 2/h^2 + exp(x_i.^2 - x_i)';

    % Constructing a sparse matrix 
    e = ones(N,1);
    L = spdiags([-e/h^2 2*e/h^2 -e/h^2],-1:1,N-1,N-1); 
    
    % Place diagonal into L along d = 0 
    d = 0; 
    L = spdiags(diag,d,L); 

    % Return the sorted eigenvalues of our discretized linear operator L 

    % Return selected option 
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

% Plots the eigenvalues with increasing N 
function plot_smallest_30_eigs(N_list, num_lowest)    
    % Create figure
    f1 = figure;
    hold on; 
    legend('FontSize',12,'Location','northwest', 'Interpreter','latex');
    
    % 30 lowest eigenvalues for N = 1,000,000 
    mill_eig = find_eigs(10^(6), num_lowest);

    % Vector of integers 
    n_vec = linspace(1,num_lowest, num_lowest);  
    plot(n_vec, mill_eig, 'o', 'Color','black','DisplayName','N=1 million')

    % Calcualtes lowest Eigenvalues in N_list 
    for i = N_list
        % Loops through and plots 30 lowest of N eigenvalues for increasing values 
        eigenvalues = find_eigs(i, num_lowest);
        plot(n_vec, eigenvalues, 'o','DisplayName', strcat('N = ', num2str(i)));
    end

    % Plot settings 
    xlabel('$n$', 'FontSize',15,'Interpreter','latex');
    ylabel('$\lambda$', 'Interpreter','latex', 'FontSize',15,'Rotation',0);
    title('Plot of 30 Smallest Eigenvalues for increasing $N$: Part 2','FontSize',18, 'Interpreter','latex');
    grid on  
end 

