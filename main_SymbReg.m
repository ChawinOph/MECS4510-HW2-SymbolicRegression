%% Assignment 2: Symbolic Regression
clc;
close all;
clear;
%% GA

down_sample_no = 10;
n_heap = 5;
n_pop = 250;
p_c = 0.8;
p_m = 0.01;
n_crossover = 2;
n_mutation = 1;
n_eval = 1e5;
n_elite = 2;

n_tour = 2;
p_tour = 0.90;

GA = GA_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
    p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour);
GA.plotXYScatter(true);
GA.plotFittest(); 
title('pre')

tic
GA.evaluate();
toc

GA.plotXYScatter(true); 
GA.plotFittest(); 
title('post')

figure
plot(GA.fittest)
