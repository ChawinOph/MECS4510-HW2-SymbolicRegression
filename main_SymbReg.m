%% Assignment 2: Symbolic Regression
clc;
close all;
clear;
%% GA

down_sample_no = 10;
n_heap = 5; % or 4?
n_pop = 250;
p_c = 0.95; 
p_m = 0.1;
n_crossover = 2;
n_mutation = 1;
n_eval = 4e5;
n_elite = 1;
trunc_rate = 1;

n_tour = 2;
p_tour = 0.90;

GP = GP_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
    p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate);
GP.plotXYScatter(true);
GP.plotFittest(); 
title('pre')

tic
GP.evaluate();
toc

GP.plotXYScatter(true); 
GP.plotFittest(); 
title('post')

figure
plot(GP.fittest)

figure 
GP.plotDot()

GP.updateFittestExpression()