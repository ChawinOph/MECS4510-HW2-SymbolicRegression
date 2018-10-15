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
        p_tour              % prob used in tournament selection
        n_elite             % number of top parents to be kept to the next generation
    end
    
    properties              % Dependent variables
        n_gen               % number of evaluations/n_pop
        fitness             % 1 x n_pop:  array of fitness correspoinding to each chromosome (mean absolute error)
        fitness_hist        % n_eval x n_pop:  array of fitness correspoinding to each chromosome over all evals
        fitness_hist_gen    % n_gen x n_pop:  array of fitness correspoinding to each chromosome over all evals
        y_func_mat          % length(this.points) x n_pop array: repeated n_pop columns of y position of points
        y_mat               % length(this.points) x n_pop array output from heaps
        fittest             % n_eval x 1 :  array of best fitness value of each eval
        fittest_gen         % n_gen x 1 :  array of best fitness value of each gen
        fittest_indv        % The current fittest individual in the pool
    end
    
    methods
        %% Constructor
        function this = GA_SymbReg(filename, down_sample_no, n_pop, n_heap, ...
                p_c, p_m, n_crossover, n_mutation, n_eval, n_tour, p_tour, n_elite, trunc_rate)
            this.filename = filename;
            this.down_sample_no = down_sample_no;
            this.points = this.importPath(this.filename);
            this.n_pop = n_pop;
            this.n_heap = n_heap;
            this.pool = this.contructHeaps();
            this.trunc_rate = trunc_rate;
            this.y_func_mat = repmat(this.points(:,2), 1, this.n_pop);
            this.updateFitness();
            
            this.p_c = p_c;
            this.p_m = p_m;
            this.n_crossover = n_crossover;
            this.n_mutation = n_mutation;
            this.n_eval = n_eval;
            this.n_gen = this.n_eval/this.n_pop;
            this.fitness_hist = zeros(this.n_eval, this.n_pop);
            this.fitness_hist_gen = zeros(this.n_gen, this.n_pop);
            this.fittest = zeros(this.n_gen, 1);
            this.fittest_gen = zeros(this.n_gen, 1);
            
            this.n_tour = n_tour;
            this.p_tour = p_tour;
            this.n_elite = n_elite;
        end
        
        %% Member Functions
        function evaluate(this)
            % evaluate: Perform an evaluation (evolving to get a new set of
            % n_pop offspring)
            for n = 1: this.n_gen
                % create a new array of offspring
                %                 offspring = nan*ones(2^this.n_heap - 1, this.n_pop);
                offspring = [];
                
                for i = 1 : ceil((this.n_pop - this.n_elite)/2)
                    % Tournament selection
                    % Randomly Select (uniform) parent candidates and osrt by fitness
                    parents = nan*ones(2^this.n_heap - 1, 2);
                    
                    for n_round = 1:2 % run through two rounds of tournaments to get two winning parents
                        cand_parent_indcs = randperm(this.n_pop, 2);
                        [~, sorted_cand_indcs] = sort(this.fitness(cand_parent_indcs));
                        if rand <= this.p_tour
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
                        
                        % find lists of heap indices that are not NaN
                        heap_indcs_parent1 = find(~isnan(parents(:, 1)));
                        heap_indcs_parent2 = find(~isnan(parents(:, 2)));
                        
                        % find limited no. of extended levels from each
                        % individual heap index in a parent
                        % index (find maximum values by running through each heap index bottom up)
                        [extended_lvl_parent1, extended_limit_parent1] = this.findExtendedLevelAndLimit(heap_indcs_parent1);
                        [extended_lvl_parent2, extended_limit_parent2] = this.findExtendedLevelAndLimit(heap_indcs_parent2);
                        
                        % Random heap index of the parent1 first (inclduing
                        % the first index for diversity?) and
                        % then check the remaining level underneath.
                        rand_idx1 = randi([1, length(heap_indcs_parent1)], 1);
                        cross_heap_indx_parent1 = heap_indcs_parent1(rand_idx1);
                        
                        % get the extesion limit from the parent
                        extended_limit_rand_parent1 = extended_limit_parent1(rand_idx1);
                        
                        % filter a set of parent2 index candidate that has the
                        % extension of the point within the limit of the
                        % random point from parent 1
                        feasible_heap_indcs_parent2 = heap_indcs_parent2(extended_lvl_parent2 <= extended_limit_rand_parent1);
                        
                        % filter a set of parent2 index candidate that has
                        % enough extension limit for the crossover point
                        % from parent1
                        feasible_heap_indcs_parent2 = feasible_heap_indcs_parent2(extended_lvl_parent1(rand_idx1) <= ...
                            extended_limit_parent2(extended_lvl_parent2 <= extended_limit_rand_parent1));
                        
                        % random a point from the feasible point list in
                        % parent 2
                        rand_idx2 = randi([1, length(feasible_heap_indcs_parent2)], 1);
                        cross_heap_indx_parent2 = feasible_heap_indcs_parent2(rand_idx2);
                        
                        parent1_extend_swapped_lvl = extended_lvl_parent1(heap_indcs_parent1 == cross_heap_indx_parent1);
                        parent2_extend_swapped_lvl = extended_lvl_parent2(heap_indcs_parent2 == cross_heap_indx_parent2);
                        parent1_extend_removed_lvl = extended_limit_parent1(heap_indcs_parent1 == cross_heap_indx_parent1);
                        parent2_extend_removed_lvl = extended_limit_parent2(heap_indcs_parent2 == cross_heap_indx_parent2);
                        
                        swapped_heap_inds_parent1 = this.findExtendedIndices(cross_heap_indx_parent1, parent1_extend_swapped_lvl);
                        swapped_heap_inds_parent2 = this.findExtendedIndices(cross_heap_indx_parent2, parent2_extend_swapped_lvl);
                        
                        replaced_heap_inds_parent1 = this.findExtendedIndices(cross_heap_indx_parent1, parent2_extend_swapped_lvl);
                        replaced_heap_inds_parent2 = this.findExtendedIndices(cross_heap_indx_parent2, parent1_extend_swapped_lvl);
                        
                        removed_heap_inds_parent1 =  this.findExtendedIndices(cross_heap_indx_parent1, parent1_extend_removed_lvl);
                        removed_heap_inds_parent2 = this.findExtendedIndices(cross_heap_indx_parent2, parent2_extend_removed_lvl);
                        
                        % get a piece of heap from parent1 that is as long
                        % as the extension in parent 2 and erase the heap
                        % all the way down to the bottom
                        swapped_heap_from_parent1 = parents(swapped_heap_inds_parent1, 1);
                        parents(removed_heap_inds_parent1, 1) = NaN;
                        % vice versa
                        swapped_heap_from_parent2 = parents(swapped_heap_inds_parent2, 2);
                        parents(removed_heap_inds_parent2, 2) = NaN;
                        
                        % put the stored swapped heaps into the other
                        % parent
                        parents(replaced_heap_inds_parent1, 1) = swapped_heap_from_parent2;
                        parents(replaced_heap_inds_parent2, 2) = swapped_heap_from_parent1;
                        
                    end
                    
                    if rand <= this.p_m
                        % mutate the parents
                        
                        % change constants by small amount
                        % swap operators
                    end
                    
                    offspring  = [offspring, parents]; %#ok<AGROW>
                end
                
                if mod(this.n_pop - this.n_elite, 2) == 1 % odd number of population
                    % randomly remove one offspring
                    remove_indx = randi(this.n_pop, 1);
                    offspring = offspring(:, 1:end ~= remove_indx);
                end
                
                % Add elites in the current generation to the current genes
                [~, sorted_fitness_indcs] = sort(this.fitness);
                offspring = [offspring, this.pool(: , sorted_fitness_indcs(1:this.n_elite))]; %#ok<AGROW>
                
                this.pool = offspring;
                this.simplify();
                fval = this.updateFitness();
                this.fitness_hist(1 + this.n_pop*(n - 1): this.n_pop*n, :) = repmat(fval, this.n_pop, 1);
                this.fitness_hist_gen(n, :) = fval;
                this.fittest(1 + this.n_pop*(n - 1): this.n_pop*n) = min(this.fitness)*ones(this.n_pop, 1);
                this.fittest_gen(n) = min(this.fitness);
                
            end
        end
        
        function extended_indcs = findExtendedIndices(~, start_index, no_extended_level)
            extended_indcs = zeros(1, 2^(no_extended_level + 1) - 1);
            extended_indcs(1) = start_index;
            for n_extend = 1 : length(extended_indcs)
                if  n_extend <= 2^(no_extended_level + 1 - 1) - 1
                    extended_indcs(2*n_extend) = 2*extended_indcs(n_extend);
                    extended_indcs(2*n_extend + 1) = 2*extended_indcs(n_extend) + 1;
                end
            end
        end
        
        function [extended_level, extend_limits] = findExtendedLevelAndLimit(this, heap_indcs_parent)
            % find maximum no. of extended levels from each
            % individual heap index in parent2
            % index (find maximum values by running through each heap index bottom up)
            extended_level = zeros(length(heap_indcs_parent), 1);
            extend_limits = zeros(length(heap_indcs_parent), 1);
            lvl_thresholds = 2.^(1 : this.n_heap) - 1;
            
            for j = length(heap_indcs_parent): -1: 1
                current_lvl = find(heap_indcs_parent(j) <= lvl_thresholds, 1, 'first');
                extend_limits(j) = this.n_heap - current_lvl;
                current_remain_lvl = zeros(length(heap_indcs_parent), 1);
                current_child_indx = heap_indcs_parent(j);
                for m = 1 : current_lvl - 1 % run until reaching the top
                    % find the list of parent current_child_lvl
                    if mod(current_child_indx, 2) == 0 % even number
                        % find parent
                        current_child_indx = (current_child_indx)/2;
                    else
                        % find parent
                        current_child_indx = (current_child_indx - 1)/2;
                    end
                    current_remain_lvl(heap_indcs_parent == current_child_indx) = m;
                end
                % get the maximum values from the current and
                % the global remain lvl
                extended_level = max([extended_level, current_remain_lvl],[],2);
            end
        end
        
        function extended_indcs = findExtendedIndicesToEnd(this, start_index)
             lvl_thresholds = 2.^(1 : this.n_heap) - 1;
             current_lvl = find(start_index <= lvl_thresholds, 1, 'first');
             extending_limits = this.n_heap - current_lvl;
             extended_indcs = this.findExtendedIndices(start_index, extending_limits);
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
        
        function m = contructHeaps(this, n_lvl)
            if nargin > 1
                m = nan*ones(2^n_lvl - 1, this.n_pop);
            else
                m = nan*ones(2^this.n_heap - 1, this.n_pop);
            end
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
        
        function simplify(this)
            % simplify (may do it every other a certain number of gens)
            % snipping: replace a sub-branch with a constant (average)
            % pruning: Eliminate sub-braches with relatively low
            % contribution (if evaluation is under a threshold,replace the heap node by a constant of zero)
            % go through each level of the heap
            min_threshold = 1e-4; % threshold for detecting minimum contribution
            denom_threshold = 1e-2; % threshold for detecting /0
            
            for i = 2^(this.n_heap) - 2 : - 2 :1
                % check if both children are both variables x or both
                % contants
                if ~isempty(find((this.pool(i, :) == 20 & this.pool(i + 1, :) == 20) | (this.pool(i, :) <= 10 & this.pool(i + 1, :) <= 10), 1))
                    indcs = find((this.pool(i, :) == 20 & this.pool(i + 1, :) == 20) | (this.pool(i, :) <= 10 & this.pool(i + 1, :) <= 10));
                    % if the parent is the minus operator
                    if ~isempty(find(this.pool(i/2, indcs) == 12, 1)) % check if there are any "-"
                        sub_indcs = indcs(this.pool(i/2, indcs) == 12);
                        if ~isempty(find(abs(this.pool(i, sub_indcs) - this.pool(i + 1, sub_indcs)) < min_threshold, 1))
                            % if there exists some columns that the have
                            % absolute values under threshold
                            subsub_indcs = sub_indcs(find(abs(this.pool(i, sub_indcs) - this.pool(i + 1, sub_indcs)) < min_threshold)); %#ok<FNDSB>
                            % change the parent to simplified values
                            this.pool(i/2, subsub_indcs) = 0;
                            % remove the children
                            this.pool(i, subsub_indcs) = NaN;
                            this.pool(i + 1, subsub_indcs) = NaN;
                        end
                    end
                    
                    % if the parent is the divide operator
                    if ~isempty(find(this.pool(i/2, indcs) == 13, 1)) % check if there are any "/"
                        sub_indcs = indcs(this.pool(i/2, indcs) == 13);
                        if ~isempty(find(abs(this.pool(i, sub_indcs )./this.pool(i + 1, sub_indcs )) == 1, 1))
                            % if there exists some columns that the have
                            % value of one because of the same values being
                            % divide
                            subsub_indcs = sub_indcs(find(this.pool(i, sub_indcs)./this.pool(i + 1, sub_indcs) == 1)); %#ok<FNDSB>
                            % change the parent to simplified values
                            this.pool(i/2, subsub_indcs) = 1;
                            % remove the children
                            this.pool(i, subsub_indcs) = NaN;
                            this.pool(i + 1, subsub_indcs) = NaN;
                        end
                    end
                    
                end
                   
                % clean up the 1 and 0
                % having parent as "/" operator
                if ~isempty(find(this.pool(i/2, :)== 13, 1))
                    % check if there exists any low values of
                    % denominator or the exact zero
                    indcs = find(this.pool(i/2, :) ==  13);
                    if ~isempty(find(abs(this.pool(i + 1, indcs)) < denom_threshold ,1))
                        % if the denominator is too small, replace by
                        % a random constant
                        sub_indcs = indcs(find(abs(this.pool(i + 1, indcs)) < denom_threshold)); %#ok<FNDSB>
                        this.pool(i + 1, sub_indcs) = 20*rand(size(sub_indcs)) - 10;
                             
                        % check if there exists any denominator = 1
                    elseif ~isempty(find(this.pool(i + 1, indcs) == 1, 1))
                        % if the denominator one, replace the "/" by the heap
                        % in the numerator heap (entire one)
                        sub_indcs = indcs(find(this.pool(i + 1, indcs) == 1)); %#ok<FNDSB>
                        
                        % overwrite heap
                        this.moveChildUp(i, sub_indcs);
%                         % find indices that will become NaN temporary
%                         removed_indcs = this.findExtendedIndicesToEnd(i/2);
%                         
%                         % find extending limit of the first child
%                         [~, extend_limit] = this.findExtendedLevelAndLimit(i);
%                        
%                         % find indcs that will be overwritten by the moving
%                         % up heap
%                         replaced_indcs = this.findExtendedIndices(i/2, extend_limit);
%                         moving_up_heap_inds = this.findExtendedIndicesToEnd(i);
%                         
%                         % store the heaps that are moving up before clearing the heaps
%                         moving_up_heaps = this.pool(moving_up_heap_inds, sub_indcs);
%                         
%                         % clear the heaps to the end
%                         this.pool(removed_indcs, sub_indcs) = NaN;
%                         
%                         % put the moving up heaps back again
%                         this.pool(replaced_indcs, sub_indcs) = moving_up_heaps;
                                                                      
                    end
                end
                
            end
        end
        
        function moveChildUp(this, first_child_indx, col_indcs)
            % Move the entire heap under i branch to i/2 branch (remove brach under i+1)
            % find indices that will become NaN temporary
            removed_indcs = this.findExtendedIndicesToEnd(first_child_indx/2);
            
            % find extending limit of the first child
            [~, extend_limit] = this.findExtendedLevelAndLimit(first_child_indx);
            
            % find indcs that will be overwritten by the moving
            % up heap
            replaced_indcs = this.findExtendedIndices(first_child_indx/2, extend_limit);
            moving_up_heap_inds = this.findExtendedIndicesToEnd(first_child_indx);
            
            % store the heaps that are moving up before clearing the heaps
            moving_up_heaps = this.pool(moving_up_heap_inds, col_indcs);
            
            % clear the heaps to the end
            this.pool(removed_indcs, col_indcs) = NaN;
            
            % put the moving up heaps back again
            this.pool(replaced_indcs, col_indcs) = moving_up_heaps;
        end
        
        function fval = updateFitness(this)
            % Evaluate the function constructed from each heap
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
            % sort individual by ascending order of MAE
            %             [this.fitness , sorted_indcs] = sort(this.fitness, 'ascend');
            %             this.pool = this.pool(:, sorted_indcs);
            [~, min_indx] = min(this.fitness);
            this.fittest_indv = this.pool(:, min_indx);
            fval = this.fitness;
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
        
        function plotDot(this)
            % Dot plot
            gens = reshape(repmat(1: this.n_gen, this.n_pop, 1), [], 1);
            scat1 = scatter(gens, reshape(this.fitness_hist_gen', [], 1));
            ylim([0 1]);
            scat1.Marker = '.';
            scat1.MarkerEdgeColor = 'b';
        end
        
    end
end