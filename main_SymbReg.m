%% Assignment 2: Symbolic Regression
clc;
close all;
clear;

%% GP_var1
down_sample_no = 10;
n_heap = 5; % or 4?
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

GP_var1 = cell(8,1);
% store best fitness value for each run
GP_var1_bestfitness = nan(length(GP_var1) ,1);
% cell array that store expressions
GP_var1_best_express = cell(length(GP_var1) ,1);
% store arrays for plotting learning curve
GP_var1_fittess_hist = nan(length(GP_var1), n_eval);

% create cell of GP
disp('GP_var1')
for i = 1:size(GP_var1_fittess_hist, 1)
    GP_var1{i} = GP_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
        p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
    tic
    GP_var1{i}.evaluate();
    toc
    GP_var1_fittess_hist(i,:) = GP_var1{i}.fittest;
    GP_var1_bestfitness(i) = GP_var1{i}.best_fitness;
    GP_var1_best_express{i} = GP_var1{i}.best_express;
end

[GP_var1_best_run_fitness_value, best_run_indx] = min(GP_var1_bestfitness);
GP_var1_best_n_eval = GP_var1{best_run_indx}.n_eval_stop;
GP_var1_best_run_express = GP_var1_best_express{best_run_indx};

figure;
GP_var1{best_run_indx}.plotXYScatter; hold on;
GP_var1{best_run_indx}.plotFittest; 

x = 0:0.1:10;
y = sin(1.5.*x)./(0.5.*x + 1);
plot(x,y)

figure;
GP_var1{best_run_indx}.plotDot; 
title('Dot Plot', 'interpreter', 'latex')
xlabel('No. of Generations','interpreter', 'latex')
ylabel('Mean Absolute Error', 'interpreter', 'latex');

avg_GP_var1 = mean(GP_var1_fittess_hist);
SEM_GP_var1 = std(GP_var1_fittess_hist)/sqrt(size(GP_var1_fittess_hist, 1));
bar_freq = n_eval/20;

save('Results\avg_GP_var1.mat','avg_GP_var1');
save('Results\SEM_GP_var1.mat','SEM_GP_var1');
save('Results\GP_var1_best_run_express.mat','GP_var1_best_run_express');
save('Results\GP_var1_best_n_eval.mat','GP_var1_best_n_eval');
save('Results\GP_var1_best_run_fitness_value.mat','GP_var1_best_run_fitness_value');

learning_fig = figure;
plotAvgSemiLogYWithErrorBar(avg_GP_var1, SEM_GP_var1, bar_freq, 'b', 1.5); hold on;
grid on; grid minor;  
title('Learning Curve', 'interpreter', 'latex')
xlabel('No. of Evaluations', 'interpreter', 'latex'); 
ylabel('Mean Absolute Error', 'interpreter', 'latex');

%% GP_var2


%% Random Search


%% 