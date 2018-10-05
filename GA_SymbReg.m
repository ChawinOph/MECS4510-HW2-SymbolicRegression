classdef GA_SymbReg < handle
    % SimpleGA Summary of this class goes here
    %   This class is designed for a TSP
    properties (Constant)
        n_bit = 1000; % chormosome length
    end
    
    properties
        filename            % file name of the given set of points
        points              % 1000 x 2 set of f(x) = x points from the text file
        down_sample_no = 1  % no. of skipped x points (use every other down_sample_no points)
        max_heap_level      % maximum number of heap levels
        n_pop               % number of population
        pool                % n_bit x n_pop array of current chormosomes
        p_c                 % cross over prob (single crossover)
        p_m                 % mutation prob
        n_crossover         % number of crossover indice of parent chromosomes (no. of cutting lines)
        n_mutation          % number of mutation indice in each parent chromosome
        n_iter              % number of evaluations/n_pop
        trunc_rate          % truncation selection ratio to improve diversity (1/n fraction form only!)
        b_find_longest      % boolean: choose goal
        n_tournament        % number of parents for competing before geeting only top two       
    end
    
    properties              % Dependent variables
        fitness             % 1 x n_pop:  array of fitness correspoinding to each chromosome
        x_mat               % n_bit x n_pop array: repeated n_pop columns of x position of points
        y_mat               % n_bit x n_pop array: repeated n_pop columns of y position of points
        fittest             % 1 x n_eval:  array of best fitness value of each iteration
        fitness_hist        % n_iter x n_pop:  array of fitness correspoinding to each chromosome over all iterations
    end
    
    methods
        %% Constuctor
        function this = GA_SymbReg(filename, down_sample_no, max_heap_size, n_pop)
            this.filename = filename;
            this.down_sample_no = down_sample_no;
            this.points = this.importPath(this.filename);
            this.max_heap_size = max_heap_size;
            this.n_pop = n_pop;
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
        
        function contructHeaps(this)
            % Randomly create the operator part
            m = zeros(this.max_heap_level, this.n_pop);
            m(1, :) = randi([1, 6], 1, this.n_pop);
            m(2 : end, :) = randi([0, 6], this.max_heap_level - 1, this.n_pop);
            
            % Change all children under the zero parents to be zeros
            
            % Reduce children under sine/cosine to be only one
            
            % Add children back to the end of each branch (if sine/cosine add only one)
            
            % Use boolean function 
            
            % Add the end of each operator as                 
        end
        
        % Snipping function
        
        
        % Pruning function
        
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