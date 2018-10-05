%% Assignment 2: Symbolic Regression
clc; 
close all;
clear;
%% GA
down_sample_no = 10;
max_heap_level = 2;
n_pop = 5;
n_crossover = 2;
n_mutation = 1;
n_eval = 1e3;

 GA = GA_SymbReg('function1.csv', down_sample_no);
 GA.plotXYScatter();

m = zeros(2^max_heap_level, n_pop);
m(1, :) = randi([1, 6], 1, n_pop);
m(2 : 2^(max_heap_level) - 2,:) = randi([0, 6], 2^(max_heap_level) - 2, n_pop);
m(end, :) = randi([0, 2], 1, n_pop);

m_logic = m == 0;




