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
        n_tour              % number of parents for competing before geeting only top (do two time)
        p_tour              
    end
    
    properties              % Dependent variables  
        n_gen               % number of evaluations/n_pop
        fitness             % 1 x n_pop:  array of fitness correspoinding to each chromosome (mean absolute error)
        fitness_hist        % n_eval x n_pop:  array of fitness correspoinding to each chromosome over all iterations
        y_func_mat          % length(this.points) x n_pop array: repeated n_pop columns of y position of points
        y_mat               % length(this.points) x n_pop array output from heaps
        fittest             % 1 x n_eval:  array of best fitness value of each iteration   
    end
    
    methods
        %% Constructor
        function this = GA_SymbReg(filename, down_sample_no, n_pop, n_heap, p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour)
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
            this.n_gen = this.n_eval/this.n_pop;
            this.fittest = min(this.fitness);
            
            this.n_tour = n_tour;
            this.p_tour = p_tour;
        end
        
        %% Member Functions
        function evaluate(this)
            % evaluate: Perform an evaluation (evolving to get a new set of
            % n_pop offspring)
            for n = 1: this.n_gen
                % create a new array of offspring
                offspring = nan*ones(2^this.n_heap - 1, this.n_pop);
                               
                for i = 1 : ceil(this.n_pop/2)
                    % Tournament selection
                    % Randomly Select (uniform) parent candidates and osrt by fitness
                    parents = nan*ones(2^this.n_heap - 1, 2);
                    for n_round = 1:2 % run through two rounds of tournaments to get two winning parents
                        cand_parent_indcs = randperm(this.n_pop, 2);
                        [~, sorted_cand_indcs] = sort(this.fitness(cand_parent_indcs));
                        if rand <= 0.9
                            % choose the fitter on
                            parents(:, n_round) = this.pool(:, cand_parent_indcs(sorted_cand_indcs(1)));
                        else
                            % choose less fit one (fluke)
                            parents(:, n_round) = this.pool(:, cand_parent_indcs(sorted_cand_indcs(2)));
                        end
                    end
                    
                    if rand <= this.p_c
                        % Crossover the parent at a randomly chosen point
                        % All children have to move along swapped parent nodes
                        
                        % check each parent of the maximum of heap index that
                        % can be chosen
                                         
                        find(~isnan(parents(:, 1)));
                        find(~isnan(parents(:, 2)));
                        
                        heap_indcs_parent1 = parents(:, 1);
                        heap_indcs_parent1(isnan(parents(:, 1))) = [];
                        heap_indcs_parent2 = parents(:, 2);
                        heap_indcs_parent2(isnan(heap_indcs_parent2)) = [];
                        
                        rand_idx1 = randperm(length(heap_indcs_parent1),1);
                        crossover_heap_indx_parent1 = heap_indcs_parent1(rand_idx1);
                        
                        rand_idx2 = randperm(length(heap_indcs_parent2),1);
                        crossover_heap_indx_parent2 = heap_indcs_parent2(rand_idx2);
                        
                        % Random heap index of the parent1 first and
                        % then check the remaining level underneath. After
                        % that, random an index of parent2 that has
                        % fewer/equal to remaining level
                        
                        % get all swapped indices of heap elements
                        cross_heap_indcs_parent1 = 2^(this.n_heap) - 2 
                        
                        cross_heap_indcs_parent2
                        
                        % prevent the heap from getting longer than the
                        % maxmimum size
                        
                        % get all swapped indices of heap elements
                        
                        swapped_heap = parents(cross_heap_indcs_parent2, 2);
                        parents(:, 2) = parents(cross_heap_indcs_parent1, 1);
                        parents(:, 1) = swapped_heap;
                        
                    end
                end
                
                
            end
        end
        
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
            this.y_mat = zeros(length(this.points), this.n_pop);
            
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
                this.y_mat(n, :) = m(1, :);
            end
            % find MAE (fitness) for each chromosome
            this.y_func_mat = repmat(this.points(:,2), 1, this.n_pop);
            this.fitness = mean(abs(this.y_func_mat - this.y_mat));
            this.fittest = mean(abs(this.y_func_mat - this.y_mat));
        end
        
        %% Visualization
        function p = plotXYScatter(this, b_hold)
            figure
            p = scatter(this.points(:,1), this.points(:,2));
            p.Marker = '.';
            p.MarkerEdgeColor = 'k';
            xlabel('X');
            ylabel('Y');
            grid on;
            grid minor;
            if nargin > 1 && b_hold
                hold on;
            end
        end
        
        function plotFittest(this)
            [~, fittest_indx] = min(this.fitness);
            plot(this.points(:,1), this.y_mat(:, fittest_indx));
        end
        
        function showSymbol(this)
            
        end
        
    end
end