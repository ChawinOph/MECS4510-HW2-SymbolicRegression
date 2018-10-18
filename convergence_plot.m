clc;
close all;
clear;

%% GP_var1
down_sample_no = 10;
n_heap = 4; % or 4?
n_pop = 250;
p_c = 0.95; 
p_m = 0.1;
n_crossover = 2;
n_mutation = 1;
n_eval = 2e5;
n_elite = 1;
trunc_rate = 1;
n_tour = 2;
p_tour = 0.90;

n_round = 4; % number of total rounds
n_run = 20;  % number of total runs per round

GP_var2 = cell(n_round, n_run);
% store arrays for plotting learning curve

percents = zeros(n_round, n_eval);

GP_var2_fittest_conv = nan(n_run, n_eval, n_round);

% create cell of GP
disp('GP_var1 convergence analysis')
for j = 1: n_round   
    for i = 1 : n_run
        GP_var2{j, i} = GP_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
            p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
        tic
        GP_var2{j, i}.evaluate();
        toc
        GP_var2_fittest_conv(i, :, j) = GP_var2{j, i}.fittest;
    end
end

convergence_MAE_threshold = 0.01;
for n = 1: n_round
    for m = 1 : n_eval
        % go through each column (20 elements from all run at a no. of evals)
        percents(n, m) = 100*sum(GP_var2_fittest_conv(:, m, n) <= convergence_MAE_threshold)/n_run;
    end
end

avg_conv_percent_GP_var2 = mean(percents);
SEM_conv_percent_GP_var2 = std(percents)/sqrt(n_round);

% save('Results\avg_conv_percent_GP_var2.mat','avg_conv_percent_GP_var2');
% save('Results\SEM_conv_percent_GP_var2.mat','SEM_conv_percent_GP_var2');
% save('Results\GP_var2_fittest_conv.mat','GP_var2_fittest_conv');

%% visualize convergence plot
bar_freq = n_eval/10;
figure;
[conv2_linehandle, conv2_e_handle] = plotAvgWithErrorBar(avg_conv_percent_GP_var2, ...
    SEM_conv_percent_GP_var2 , bar_freq, 'b', 1.5); hold on
[conv1_linehandle, conv1_e_handle] = plotAvgWithErrorBar(avg_conv_percent_GP_var1, ...
    SEM_conv_percent_GP_var1 , bar_freq, 'r', 1.5); hold on
grid on; grid minor;  
title('Convergence Plot', 'interpreter', 'latex');
xlabel('No. of Evaluations', 'interpreter', 'latex'); 
ylabel(['$\%$ of Runs with MAE $\le$ ', num2str(convergence_MAE_threshold)], 'interpreter', 'latex');

legnd = legend([conv2_linehandle, conv1_linehandle, conv2_e_handle, conv1_e_handle], ...
   {'Genetics Programming Max. Heap Level = 4',...
    'Genetics Programming Max. Heap Level = 5' ,...
    '$\mu\pm1\sigma_{\bar{x}}$','$\mu\pm1\sigma_{\bar{x}}$'});
set(legnd,'Interpreter','latex')
legnd.NumColumns = 2;
