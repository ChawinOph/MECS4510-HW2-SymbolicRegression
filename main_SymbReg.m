%% Assignment 2: Symbolic Regression
clc;
close all;
clear;
%% GA

down_sample_no = 10;
n_heap = 3;
n_pop = 10;
p_c = 0.5;
p_m = 0.01;
n_crossover = 2;
n_mutation = 1;
n_eval = 1e3;

n_tour = 2;
p_tour = 0.9;

GA = GA_SymbReg('function1.csv',down_sample_no, n_pop, n_heap,...
    p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour);
GA.plotXYScatter(true);
GA.plotFittest(); 
GA.evaluate();

%%

m = nan*ones(2^n_heap - 1, n_pop);
m(1, :) = randi([11, 16], 1, n_pop);
for i = 1 : 2^(n_heap) - 1
    % Running through each index of the heap structure
    if ~isempty(find(m(i, :)  == 11 | m(i, :) == 14, 1))
        % operators +, -, /, *
        if i <= 2^(n_heap - 2) - 1
            m([2*i, 2*i + 1], m(i, :)  == 11 | m(i, :) == 14) = randi([11, 18], ...
                2, size(m(m(i, :)  == 11 | m(i, :) == 14), 2));
        else
            % the row before the last row can have x or constant children
            m([2*i, 2*i + 1], m(i, :)  == 11 | m(i, :) == 14) = randi([17, 18], ...
                2, size(m(m(i, :)  == 11 | m(i, :) == 14), 2));
        end
    end
    
    if ~isempty(find(m(i, :)  == 12 | m(i, :) == 13, 1))
        % operators +, -, /, *
        if i <= 2^(n_heap - 2) - 1
            m([2*i, 2*i + 1], m(i, :)  == 12 | m(i, :) == 13) = randi([11, 18], ...
                2, size(m(m(i, :)  == 12 | m(i, :) == 13), 2));
        else
            % the row before the last row can have x or constant children
            m([2*i, 2*i + 1], m(i, :)  == 12 | m(i, :) == 13) = randi([17, 18], ...
                2, size(m(m(i, :)  == 12 | m(i, :) == 13), 2));
        end
    end
    
    if ~isempty(find(m(i, :) >= 15 & m(i, :) <= 16, 1))
        % sin or cosine cases
        if i <= 2^(n_heap - 2) - 1 % not the row before the last row
            m(2*i, m(i, :) >= 15 & m(i, :) <= 16) = randi([11, 18], ...
                1, size(m(m(i, :) >= 15 & m(i, :) <= 16), 2));
        elseif i <= 2^(n_heap - 1) - 1 % the row before the last row can only have x or constant children
            m(2*i, m(i, :) >= 15 & m(i, :) <= 16) = randi([17, 18], ...
                1, size(m(m(i, :) >= 15 & m(i, :) <= 16), 2));
        end
        m(2*i + 1, m(i, :) >= 15 & m(i, :) <= 16) = nan;
    end
    
    if ~isempty(find(m(i, :) == 17, 1))
        % If the index is a terminal of a constant
        % assign all zeros underneath an assign a contant at the 7
        m(i, m(i, :) == 17) = 20*rand(size(find(m(i, :) == 17))) - 10;
        if i <= 2^(n_heap - 2) - 1
            % if still not the row before the last row, assign zeros
            % underneath
            m([2*i, 2*i + 1], m(i, :) == 17) = nan;
        end
    end
    
    if ~isempty(find(m(i, :) == 18, 1))
        % If the index is a terminal (constant or x)
        % assign all zeros underneath an assign a contant at the 7
        m(i, m(i, :) == 18) = 20; % change to a number beyond the constant range
        if i <= 2^(n_heap - 2) - 1
            % if still not the row before the last row, assign zeros
            % underneath
            %             m([2*i, 2*i + 1], m(i, :) == 18) = zeros(2, size(m(m(i, :) == 18), 2));
            m([2*i, 2*i + 1], m(i, :) == 18) = nan;
        end
    end
    
end

GA.pool = m;

%% evaluate the matrix
% replace 20's by a value of x
y_mat = zeros(length(GA.points), GA.n_pop);

for n = 1: length(GA.points)
    m = GA.pool;
    m(m == 20) = GA.points(n, 1);
    
    for i = 2^(n_heap) - 2 : - 2 :1
        if ~isempty(find(m(i/2, :) == 11, 1)) % check if there are any "+"
            m(i/2 , m(i/2, :) == 11) = m(i, m(i/2, :) == 11) + m(i + 1, m(i/2, :) == 11);
        end
        
        if ~isempty(find(m(i/2, :) == 12, 1)) % check if there are any "-"
            m(i/2 , m(i/2, :) == 12) = m(i, m(i/2, :) == 12) - m(i + 1, m(i/2, :) == 12);
        end
        
        if ~isempty(find(m(i/2, :) == 13, 1)) % check if there are any "/"
            m(i/2 , m(i/2, :) == 13) = m(i, m(i/2, :) == 13)./m(i + 1, m(i/2, :) == 13);
        end
        
        if ~isempty(find(m(i/2, :) == 14, 1)) % check if there are any "*"
            m(i/2 , m(i/2, :) == 14) = m(i, m(i/2, :) == 14).*m(i + 1, m(i/2, :) == 14);
        end
        
        if ~isempty(find(m(i/2, :) == 15, 1)) % check if there are any "sine"
            m(i/2 , m(i/2, :) == 15) = sin(m(i, m(i/2, :) == 15));
        end
        
        if ~isempty(find(m(i/2, :) == 16, 1)) % check if there are any "cosine"
            m(i/2 , m(i/2, :) == 16, :) = cos(m(i, m(i/2, :) == 16));
            m([i, i + 1], m(i/2, :) == 16) = nan;
        end
    end
    y_mat(n, :) = m(1, :);
end

% find MAE (fitness) for each chromosome
y_func_mat = repmat(GA.points(:,2), 1, GA.n_pop);
mean_abs_err = mean(abs(y_func_mat - y_mat));