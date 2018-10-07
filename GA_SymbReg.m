classdef GA_SymbReg < handle
    % SimpleGA Summary of this class goes here
    %   This class is designed for a TSP
    properties
        filename            % file name of the given set of points
        points              % 1000 x 2 set of f(x) = x points from the text file
        down_sample_no      % no. of skipped x points (use every other down_sample_no points)
        n_heap              % maximum number of heap levels
        n_pop               % number of population
        pool                % n_heap  x n_pop array of current chormosomes (heaps)
        
        p_c                 % cross over prob (single crossover)
        p_m                 % mutation prob
        n_crossover         % number of crossover indice of parent chromosomes (no. of cutting lines)
        n_mutation          % number of mutation indice in each parent chromosome
        n_eval              % number of evaluations
        
        trunc_rate          % truncation selection ratio to improve diversity (1/n fraction form only!
        n_tournament        % number of parents for competing before geeting only top two
    end
    
    properties              % Dependent variables  
        n_gen               % number of evaluations/n_pop
        fitness             % 1 x n_pop:  array of fitness correspoinding to each chromosome (mean absolute error)
        fitness_hist        % n_eval x n_pop:  array of fitness correspoinding to each chromosome over all iterations
        y_func_mat          % length(this.points) x n_pop array: repeated n_pop columns of y position of points
        fittest             % 1 x n_eval:  array of best fitness value of each iteration   
    end
    
    methods
        %% Constructor
        function this = GA_SymbReg(filename, down_sample_no, n_pop, n_heap, p_c, p_m, n_crossover, n_mutation, n_eval, n_gen)
            this.filename = filename;
            this.down_sample_no = down_sample_no;
            this.points = this.importPath(this.filename);
            this.n_pop = n_pop;
            this.n_heap = n_heap;
            this.pool = this.contructHeaps();
            this.y_func_mat = repmat(this.points(:,2), 1, this.n_pop);
            this.updateFitness();
            
            this.p_c = p_c;
            this.p_m = p_m;
            this.n_crossover = n_crossover;
            this.n_mutation = n_mutation;
            this.n_eval = n_eval;
            this.n_gen = n_gen;
        end
        
        %% Member Functions
        function points = importPath(this, csv_filename)
            % Import the text file or array
            if ischar(csv_filename)
                dlist = dir(csv_filename);
                file_directory = [dlist.folder '\' dlist.name] ;
                points = csvread(file_directory);
            elseif isfloat(csv_filename)
                points = csv_filename;
            end
            points = points(1: this.down_sample_no: end, :);
        end
        
        function m = contructHeaps(this)
            m = nan*ones(2^this.n_heap - 1, this.n_pop);
            m(1, :) = randi([11, 16], 1, this.n_pop);
            for i = 1 : 2^(this.n_heap) - 1
                % Running through each index of the heap structure
                if ~isempty(find(m(i, :)  == 11 | m(i, :) == 14, 1))
                    % operators +, -, /, *
                    if i <= 2^(this.n_heap - 2) - 1
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
                    if i <= 2^(this.n_heap - 2) - 1
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
                    if i <= 2^(this.n_heap - 2) - 1 % not the row before the last row
                        m(2*i, m(i, :) >= 15 & m(i, :) <= 16) = randi([11, 18], ...
                            1, size(m(m(i, :) >= 15 & m(i, :) <= 16), 2));
                    elseif i <= 2^(this.n_heap - 1) - 1 % the row before the last row can only have x or constant children
                        m(2*i, m(i, :) >= 15 & m(i, :) <= 16) = randi([17, 18], ...
                            1, size(m(m(i, :) >= 15 & m(i, :) <= 16), 2));
                    end
                    m(2*i + 1, m(i, :) >= 15 & m(i, :) <= 16) = nan;
                end
                
                if ~isempty(find(m(i, :) == 17, 1))
                    % If the index is a terminal of a constant
                    % assign all zeros underneath an assign a contant at the 7
                    m(i, m(i, :) == 17) = 20*rand(size(find(m(i, :) == 17))) - 10;
                    if i <= 2^(this.n_heap - 2) - 1
                        % if still not the row before the last row, assign zeros
                        % underneath
                        m([2*i, 2*i + 1], m(i, :) == 17) = nan;
                    end
                end
                
                if ~isempty(find(m(i, :) == 18, 1))
                    % If the index is a terminal (constant or x)
                    % assign all zeros underneath an assign a contant at the 7
                    m(i, m(i, :) == 18) = 20; % change to a number beyond the constant range
                    if i <= 2^(this.n_heap - 2) - 1
                        % if still not the row before the last row, assign zeros
                        % underneath
                        %             m([2*i, 2*i + 1], m(i, :) == 18) = zeros(2, size(m(m(i, :) == 18), 2));
                        m([2*i, 2*i + 1], m(i, :) == 18) = nan;
                    end
                end
            end
        end
        
        function updateFitness(this)
            y_mat = zeros(length(this.points), this.n_pop);
            
            for n = 1: length(this.points)
                m = this.pool;
                m(m == 20) = this.points(n, 1);
                
                for i = 2^(this.n_heap) - 2 : - 2 :1
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
            this.y_func_mat = repmat(this.points(:,2), 1, this.n_pop);
            this.fitness = mean(abs(this.y_func_mat - y_mat));
        end
        
        %% Visualization
        function p = plotXYScatter(this)
            figure
            p = scatter(this.points(:,1), this.points(:,2));
            p.Marker = '.';
            p.MarkerEdgeColor = 'k';
            xlabel('X');
            ylabel('Y');
            grid on;
            grid minor;
        end
        
        function showSymbol(this)
            
        end
        
    end
end