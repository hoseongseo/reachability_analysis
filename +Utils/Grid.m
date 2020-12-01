classdef Grid < handle
    %%%% Unit grid points on n dimensional space
    %%%% x_i \in [0,1] for all i = 1:n
    properties (GetAccess = public, SetAccess = public)
        x       % points
    end
    methods (Static)
        function this = Grid(numPts)
            % generates points in [0,1]
            % numPts : number of points in each dimension
            n = length(numPts);
            
            this.x = cell(1,n);
            gr = cell(1,n);
            for i = 1:n
                gr{i} = linspace(0,1,numPts(i));
            end
            
            if n == 1
                this.x{1} = gr{1}';
            else
                for i = 1:n
                    n_gr = length(gr{i});
                    
                    sz_basis = ones(1,n);
                    sz_basis(i) = n_gr;
                    
                    sz_all = numPts;
                    sz_all(i) = 1;
                    
                    basis = reshape(gr{i}, sz_basis);
                    this.x{i} = repmat(basis, sz_all);
                end
            end
        end
        function this = rand(numPts)
            % generates random points in [0,1]
            % numPts : number of points in each dimension
            n = length(numPts);
            
            this.x = cell(1,n);
            
            if n == 1
                this.x{1} = rand(numPts,1);
            else
                for i = 1:n
                    this.x{i} = rand( numPts );
                end
            end
        end
    end
end