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

GP_var1 = cell(5,1);
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
title('Dataset (down sampled) and Best Fit Curve', 'interpreter', 'latex')
leg = legend('Dataset', 'Best Fit Curve'); set(leg,'Interpreter','latex')

% x = 0:0.1:10;
% y = sin(1.5.*x)./(0.5.*x + 1);
% plot(x,y)

figure;
GP_var1{best_run_indx}.plotDot; 
title('Dot Plot', 'interpreter', 'latex')
xlabel('No. of Generations','interpreter', 'latex')
ylabel('Mean Absolute Error', 'interpreter', 'latex');

avg_GP_var1 = mean(GP_var1_fittess_hist);
SEM_GP_var1 = std(GP_var1_fittess_hist)/sqrt(size(GP_var1_fittess_hist, 1));
bar_freq = n_eval/10;

% save('Results\avg_GP_var1.mat','avg_GP_var1');
% save('Results\SEM_GP_var1.mat','SEM_GP_var1');
% save('Results\GP_var1_best_run_express.mat','GP_var1_best_run_express');
% save('Results\GP_var1_best_n_eval.mat','GP_var1_best_n_eval');
% save('Results\GP_var1_best_run_fitness_value.mat','GP_var1_best_run_fitness_value');

%% GP_var2 n_heap = 4;
n_heap = 4; % or 4?
GP_var2 = cell(5,1);
% store best fitness value for each run
GP_var2_bestfitness = nan(length(GP_var2) ,1);
% cell array that store expressions
GP_var2_best_express = cell(length(GP_var2) ,1);
% store arrays for plotting learning curve
GP_var2_fittess_hist = nan(length(GP_var2), n_eval);

% create cell of GP
disp('GP_var2')
for i = 1:size(GP_var2_fittess_hist, 1)
    GP_var2{i} = GP_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
        p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
    tic
    GP_var2{i}.evaluate();
    toc
    GP_var2_fittess_hist(i,:) = GP_var2{i}.fittest;
    GP_var2_bestfitness(i) = GP_var2{i}.best_fitness;
    GP_var2_best_express{i} = GP_var2{i}.best_express;
end

avg_GP_var2 = mean(GP_var2_fittess_hist);
SEM_GP_var2 = std(GP_var2_fittess_hist)/sqrt(size(GP_var2_fittess_hist, 1));

% save('Results\avg_GP_var2.mat','avg_GP_var2');
% save('Results\SEM_GP_var2.mat','SEM_GP_var2');
%% Random Search
GP_RS = cell(5,1);
% store best fitness value for each run
GP_RS_bestfitness = nan(length(GP_RS) ,1);
% cell array that store expressions
GP_RS_best_express = cell(length(GP_RS) ,1);
% store arrays for plotting learning curve
GP_RS_fittess_hist = nan(length(GP_RS), n_eval);

n_heap = 4; % or 4?

disp('Random Search')
for i = 1:size(GP_RS_fittess_hist, 1)
    GP_RS{i} = GP_SymbReg('function1.csv',down_sample_no, 1, n_heap,...
        p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
    tic
    GP_RS{i}.random_search();
    toc
    GP_RS_fittess_hist(i,:) = GP_RS{i}.fittest_gen;
end
    
avg_GP_RS = mean(GP_RS_fittess_hist);
SEM_GP_RS = std(GP_RS_fittess_hist)/sqrt(size(GP_RS_fittess_hist, 1));

% save('Results\avg_GP_RS.mat','avg_GP_RS');
% save('Results\SEM_GP_RS.mat','SEM_GP_RS');

%% Random Mutation Hill Climber (RMHC)
GP_RMHC = cell(5,1);
% store best fitness value for each run
GP_RMHC_bestfitness = nan(length(GP_RMHC) ,1);
% cell array that store expressions
GP_RMHC_best_express = cell(length(GP_RMHC) ,1);
% store arrays for plotting learning curve
GP_RMHC_fittess_hist = nan(length(GP_RMHC), n_eval);

n_heap = 4; 

disp('Random Mutation Hill Climber')
for i = 1:size(GP_RS_fittess_hist, 1)
    GP_RMHC{i} = GP_SymbReg('function1.csv',down_sample_no, 1, n_heap,...
        p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
    tic
    GP_RMHC{i}.random_mutation_hill_climber();
    toc
    GP_RMHC_fittess_hist(i,:) = GP_RMHC{i}.fittest_gen;
end
    
avg_GP_RMHC = mean(GP_RMHC_fittess_hist);
SEM_GP_RMHC = std(GP_RMHC_fittess_hist)/sqrt(size(GP_RMHC_fittess_hist, 1));

% save('Results\avg_GP_RMHC.mat','avg_GP_RMHC');
% save('Results\SEM_GP_RMHC.mat','SEM_GP_RMHC');

%% Plot Learning Curves
figure;
[p_RS, f_RS] = plotAvgSemiLogYWithErrorBar(avg_GP_RS, SEM_GP_RS, bar_freq, 'k', 1.5); hold on;
[p_RMHC, f_RMHC] = plotAvgSemiLogYWithErrorBar(avg_GP_RMHC, SEM_GP_RMHC, bar_freq, 'm', 1.5); hold on;
[p_GAvar2, f_GAvar2] = plotAvgSemiLogYWithErrorBar(avg_GP_var2, SEM_GP_var2, bar_freq, 'b', 1.5); hold on;
[p_GAvar1, f_GAvar1] = plotAvgSemiLogYWithErrorBar(avg_GP_var1, SEM_GP_var1, bar_freq, 'r', 1.5); hold on;

grid on; grid minor;  
title('Learning Curves', 'interpreter', 'latex');
xlabel('No. of Evaluations', 'interpreter', 'latex'); 
ylabel('Mean Absolute Error (MAE)', 'interpreter', 'latex');
legnd = legend([p_RS, p_RMHC, p_GAvar2, p_GAvar1, f_RS, f_RMHC, f_GAvar2, f_GAvar1],... 
   {'Random Search','Random Mutation Hill Climber','Genetics Programming Max. Heap Level = 4',...
    'Genetics Programming Max. Heap Level = 5' ,...
    '$\pm1\sigma_{\bar{x}}$','$\pm1\sigma_{\bar{x}}$','$\pm1\sigma_{\bar{x}}$', ...
    '$\pm1\sigma_{\bar{x}}$'});
set(legnd,'Interpreter','latex')
legnd.NumColumns = 2;

%% simpler problems tested for debugging
% figure;

% simple cosine function
x1 = 0:0.05:10;
y1 = cos(2.*x1);

% simple linear function
x2 = 0:0.1:10;
y2 = 2.*x2 + 1;

% simple parabolic function
% x3 = 0:0.1:10;
% y3 = -3*x3.*x3;

down_sample_no = 1;
n_heap = 4;
n_pop = 100;
p_c = 0.95;
p_m = 0.1;
n_crossover = 2;
n_mutation = 1;
n_eval = 1e4;
n_elite = 1;
trunc_rate = 1;

GP_var1_simple_cos = GP_SymbReg([x1',y1'],down_sample_no, n_pop, n_heap,...
    p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);

GP_var1_simple_linear = GP_SymbReg([x2',y2'],down_sample_no, n_pop, n_heap,...
    p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);

% GP_var1_simple_para = GP_SymbReg([x3',y3'],down_sample_no, n_pop, n_heap,...
%     p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);

[~, indx] = min(GP_var1_simple_para.fitness);

tic
GP_var1_simple_cos.evaluate()
toc

tic
GP_var1_simple_linear.evaluate()
toc

tic
GP_var1_simple_para.evaluate()
toc

figure;
% subplot(1,2,1);

GP_var1_simple_cos().plotXYScatter; hold on;
GP_var1_simple_cos().plotFittest;
title('Test Function: $f_1(x) = \cos(2x)$', 'interpreter', 'latex')
leg = legend('Data Set','Best Fit Curve');
set(leg,'Interpreter','latex')

% subplot(1,2,2);
figure;
GP_var1_simple_linear().plotXYScatter; hold on;
GP_var1_simple_linear().plotFittest;
title('Test Function: $f_2(x) = 2x + 1$', 'interpreter', 'latex')

% subplot(1,3,3);
% 
% GP_var1_simple_para().plotXYScatter; hold on;
% GP_var1_simple_para().plotFittest;
% title('Test Function: $f_3(x) = -3x^2$', 'interpreter', 'latex')

figure;
plot(GP_var1_simple_cos.fittest,'Linewidth', 1.5); hold on;
plot(GP_var1_simple_linear.fittest,'Linewidth', 1.5);
% plot(GP_var1_simple_para.fittest,'Linewidth', 1.5);
grid on; grid minor;  
title('Learning Curves using GP', 'interpreter', 'latex');
xlabel('No. of Evaluations', 'interpreter', 'latex'); 
ylabel('Mean Absolute Error (MAE)', 'interpreter', 'latex');
legnd = legend({'$f_1(x) = \cos(2x)$', '$f_2(x) = 2x + 1$'});
set(legnd,'Interpreter','latex')
set(gca,'YScale','log');

